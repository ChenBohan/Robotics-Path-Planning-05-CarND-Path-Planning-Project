#ifndef __VEHICLE__
#define __VEHICLE__

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Help.h"
#include "Types.h"
#include "Constants.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

class Vehicle {
    public:
        Vehicle(vector<double> map_x, vector<double> map_y, 
                vector<double> map_s, vector<double> map_dx, 
                vector<double> map_dy){
                map_waypoints_s = map_s;
                map_waypoints_x = map_x;
                map_waypoints_y = map_y;
                map_waypoints_dx = map_dx;
                map_waypoints_dy = map_dy;

                num_waypoints = map_waypoints_x.size();

                too_close = false;
            }
        void Update(CarLocalizationData& vLocal, 
            vector<double> previous_path_x, vector<double> previous_path_y,
            vector<vector<double>> sensorData){
                coarse_waypoints_s.clear();
                coarse_waypoints_x.clear();
                coarse_waypoints_y.clear();
                coarse_waypoints_dx.clear();
                coarse_waypoints_dy.clear();
                interpolated_waypoints_s.clear();
                interpolated_waypoints_x.clear();
                interpolated_waypoints_y.clear();
                interpolated_waypoints_dx.clear();
                interpolated_waypoints_dy.clear();
                sensor_fusion.clear();
                car_x = vLocal.x;
                car_y = vLocal.y;
                car_s = vLocal.s;
                car_d = vLocal.d;
                car_yaw = vLocal.yaw;
                car_speed = vLocal.speed * 0.44704;
                sensor_fusion = sensorData;

                int next_waypoint_index = NextWaypoint(car_x, car_y, car_yaw, map_waypoints_x, map_waypoints_y);
                
                for (int i = -NUM_WAYPOINTS_BEHIND; i < NUM_WAYPOINTS_AHEAD; i++) {
                    // for smooting, take so many previous and so many subsequent waypoints
                    int idx = (next_waypoint_index+i) % num_waypoints;
                    if (idx < 0) {
                        // correct for wrap
                        idx += num_waypoints;
                    }
                    // correct for wrap in s for spline interpolation (must be continuous)
                    double current_s = map_waypoints_s[idx];
                    double base_s = map_waypoints_s[next_waypoint_index];
                    if (i < 0 && current_s > base_s) {
                        current_s -= TRACK_LENGTH;
                    }
                    if (i > 0 && current_s < base_s) {
                        current_s += TRACK_LENGTH;
                    }
                    coarse_waypoints_s.push_back(current_s);
                    coarse_waypoints_x.push_back(map_waypoints_x[idx]);
                    coarse_waypoints_y.push_back(map_waypoints_y[idx]);
                    coarse_waypoints_dx.push_back(map_waypoints_dx[idx]);
                    coarse_waypoints_dy.push_back(map_waypoints_dy[idx]);
                }

                double dist_inc = 0.5;	
                int num_interpolation_points = (coarse_waypoints_s[coarse_waypoints_s.size()-1] - coarse_waypoints_s[0]) / dist_inc;
                // interpolated s is simply...
                interpolated_waypoints_s.push_back(coarse_waypoints_s[0]);
                for (int i = 1; i < num_interpolation_points; i++) {
                    interpolated_waypoints_s.push_back(coarse_waypoints_s[0] + i * dist_inc);
                }
                interpolated_waypoints_x = interpolate_points(coarse_waypoints_s, coarse_waypoints_x, dist_inc, num_interpolation_points);
                interpolated_waypoints_y = interpolate_points(coarse_waypoints_s, coarse_waypoints_y, dist_inc, num_interpolation_points);
                interpolated_waypoints_dx = interpolate_points(coarse_waypoints_s, coarse_waypoints_dx, dist_inc, num_interpolation_points);
                interpolated_waypoints_dy = interpolate_points(coarse_waypoints_s, coarse_waypoints_dy, dist_inc, num_interpolation_points);

                subpath_size = min(PREVIOUS_PATH_POINTS_TO_KEEP, (int)previous_path_x.size());
                traj_start_time = subpath_size * PATH_DT;
                double pos_x, pos_y, pos_x2, pos_y2, angle, vel_x1, vel_y1,
                    pos_x3, pos_y3, vel_x2, vel_y2, acc_x, acc_y;

                if (subpath_size < 4) {
                    pos_x = car_x;
                    pos_y = car_y;
                    angle = deg2rad(car_yaw);
                    pos_s = car_s;
                    pos_d = car_d;
                    s_dot = car_speed;
                    d_dot = 0;
                    s_ddot = 0;
                    d_ddot = 0;
                } else {
                    // consider current position to be last point of previous path to be kept
                    pos_x = previous_path_x[subpath_size-1];
                    pos_y = previous_path_y[subpath_size-1];
                    pos_x2 = previous_path_x[subpath_size-2];
                    pos_y2 = previous_path_y[subpath_size-2];
                    angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
                    vector<double> frenet = getFrenet(pos_x, pos_y, angle, interpolated_waypoints_x, interpolated_waypoints_y, interpolated_waypoints_s);
                    pos_s = frenet[0];
                    pos_d = frenet[1];

                    int next_interp_waypoint_index = NextWaypoint(pos_x, pos_y, angle, interpolated_waypoints_x, 
                                                                                                                interpolated_waypoints_y);
                    double dx = interpolated_waypoints_dx[next_interp_waypoint_index - 1];
                    double dy = interpolated_waypoints_dy[next_interp_waypoint_index - 1];
                    // sx,sy vector is perpendicular to dx,dy
                    double sx = -dy;
                    double sy = dx;

                    // calculate s_dot & d_dot
                    vel_x1 = (pos_x - pos_x2) / PATH_DT;
                    vel_y1 = (pos_y - pos_y2) / PATH_DT;
                    s_dot = vel_x1 * sx + vel_y1 * sy;
                    d_dot = vel_x1 * dx + vel_y1 * dy;

                    // have to get another point to calculate s_ddot, d_ddot from xy acceleration
                    pos_x3 = previous_path_x[subpath_size-3];
                    pos_y3 = previous_path_y[subpath_size-3];
                    vel_x2 = (pos_x2 - pos_x3) / PATH_DT;
                    vel_y2 = (pos_y2 - pos_y3) / PATH_DT;
                    acc_x = (vel_x1 - vel_x2) / PATH_DT;
                    acc_y = (vel_y1 - vel_y2) / PATH_DT;
                    s_ddot = acc_x * sx + acc_y * sy;
                    d_ddot = acc_x * dx + acc_y * dy;		
                }	
                genPredictionsFromSensorData();
                genTraj();
            }
        void testDisplay(){
                cout 
                // << " coarse_waypoints_s:" << coarse_waypoints_s.size()
                // << " coarse_waypoints_x:" << coarse_waypoints_x.size()
                // << " coarse_waypoints_y:" << coarse_waypoints_y.size()
                // << " inter_wapoints_s:" << interpolated_waypoints_s.size()
                // << " inter_wapoints_x:" << interpolated_waypoints_x.size()
                // << " inter_wapoints_y:" << interpolated_waypoints_y.size()
                // << " inter_wapoints_dx:" << interpolated_waypoints_dx.size()
                // << " inter_wapoints_dy:" << interpolated_waypoints_dy.size()
                // << " s:" << pos_s << " s_d:" << s_dot << "s_dd:" << s_ddot
                // << " d:" << pos_d << " d_d:" << s_dot << "d_dd:" << s_ddot
                // << " otherCar size:" << otherCarPredictions.size()
                // << " otherCar sd size:" << otherCarPredictions[0].size()
                << " best_state:" << best_traj_state
                << " best_cost:" << best_cost
                // << " xval_size:" << next_x_vals.size()
                // << " yval_size:" << next_y_vals.size()
                << " too close:" << too_close
                << " foolow speed:" << follow_speed
                << endl;
            }
        void genPredictionsFromSensorData(){
                // ********************* GENERATE PREDICTIONS FROM SENSOR FUSION DATA **************************
                // The data format for each car is: [ id, x, y, vx, vy, s, d]. The id is a unique identifier for that car. The x, y values are in global map coordinates, and the vx, vy values are the velocity components, also in reference to the global map. Finally s and d are the Frenet coordinates for that car.
                otherCarPredictions.clear();
                carsInLane0.clear();
                carsInLane1.clear();
                carsInLane2.clear();
                too_close = false;
                duration = N_SAMPLES * DT - subpath_size * PATH_DT;
                for (auto sf: sensor_fusion) {
                    double other_car_vel = sqrt(pow((double)sf[3], 2) + pow((double)sf[4], 2));
                    double v_s = sf[5];
                    double v_d = sf[6];
                    int v_lane = DToLane(v_d);
                    if(v_lane==0) carsInLane0.push_back(sf);
                    if(v_lane==1) carsInLane1.push_back(sf);
                    if(v_lane==2) carsInLane2.push_back(sf);
                    vector<vector<double>> predictions;
                    for (int i = 0; i < N_SAMPLES; i++){
                        double t = traj_start_time + (i * duration / N_SAMPLES);
                        double new_s = v_s + other_car_vel * t;
                        vector<double> s_and_d = {new_s, v_d};
                        predictions.push_back(s_and_d);
                    }
                    otherCarPredictions.push_back(predictions);
                }
                cout << "lane0:" << carsInLane0.size()
                    << " lane1:" << carsInLane1.size()
                    << " lane2:" << carsInLane2.size()
                    << endl;
                checkTooClose();
            }
        void checkTooClose(){
                int my_lane = DToLane(car_d);
                if(my_lane==0){
                    too_close = checkAhead(carsInLane0);
                }else if(my_lane==1){
                    too_close = checkAhead(carsInLane1);
                }else if(my_lane==2){
                    too_close = checkAhead(carsInLane2);
                }
            }
        bool checkAhead(vector<vector<double>> cars){
                double gap = 99999; 
                for(int i=0; i<cars.size(); i++){
                    double check_car_s = cars[i][5];
                    double check_car_vx = cars[i][3];
                    double check_car_vy = cars[i][4];
                    if(check_car_s>car_s){
                        double temp = check_car_s - car_s;
                        if(temp < gap){
                            gap = temp;
                            follow_speed = sqrt(pow(check_car_vx, 2) + pow(check_car_vy, 2));                            
                        }
                    }
                }
                if(gap < 15){
                    return true;
                }
                else{
                    return false;
                }
            }
        void laneChange(){
                car_to_left = false;
                car_to_right = false;
                car_just_ahead = false;
                for (vector<vector<double>> other_car: otherCarPredictions) {
                    if(other_car[0][0] > car_s){
                        double s_diff = other_car[0][0] - car_s;
                        if (s_diff < FOLLOW_DISTANCE) {
                            // cout << "s diff: " << s_diff << endl;
                            double d_diff = other_car[0][1] - car_d;
                            if (d_diff > 2 && d_diff < 6) {
                                car_to_right = true;
                            } else if (d_diff < -2 && d_diff > -6) {
                                car_to_left = true;
                            } else if (d_diff > -2 && d_diff < 2) {
                                car_just_ahead = true;
                            }
                        }
                    }
                    
                }
            }
        void updateState(){
                available_states.clear();
                available_states = {"KL"};
                if (car_d > 4 && !car_to_left) {
                    available_states.push_back("LCL");
                }
                if (car_d < 8 && !car_to_right) {
                    available_states.push_back("LCR");
                }
            }
        void genTraj(){
                laneChange();            
                updateState();
                // cout << "available_states:" << available_states.size() << endl;
                best_cost = 999999;
                best_traj_state = "";
                for (string state: available_states) {
                    vector<vector<double>> target_s_and_d = genTargetSD(state);
                    genBestTarget(target_s_and_d, state);
                }
            }
        vector<vector<double>> genTargetSD(string state) {
                double target_d, target_d_d, target_d_dd;
                double target_s, target_s_d, target_s_dd;
                current_lane = car_d / 4;
                target_d = 0;
                target_d_d = 0;
                target_d_dd = 0;
                target_s_d = min(s_dot + MAX_INSTANTANEOUS_ACCEL/4 * duration, SPEED_LIMIT);
                // target_s_d = SPEED_LIMIT;  
                target_s_dd = 0;
                target_s = car_s + (s_dot + target_s_d) / 2 * duration;

                vector<double> leading_vehicle_s_and_sdot;

                if(state.compare("KL") == 0){
                    target_d = (double)current_lane * 4 + 2;
                    target_lane = target_d / 4;
                }else if(state.compare("LCL") == 0){
                    target_d = ((double)current_lane - 1) * 4 + 2;
                    target_lane = target_d / 4;
                }else if(state.compare("LCR") == 0){
                    target_d = ((double)current_lane + 1) * 4 + 2;
                    target_lane = target_d / 4;}
                if (car_just_ahead) {
                    target_s_d = 0.0;
                }
                return {{target_s, target_s_d, target_s_dd}, {target_d, target_d_d, target_d_dd}};
            }
        void genBestTarget(vector<vector<double>> target, string state){
                vector<double> s_traj;
                vector<double> d_traj;
                vector<double> target_s = target[0];
                vector<double> target_d = target[1];
                vector<double> current_s = {this->car_s, this->s_dot, this->s_ddot};
                vector<double> current_d = {this->car_d, this->d_dot, this->d_ddot};
                this->s_traj_coeffs = calJMT(current_s, target_s, duration);
                this->d_traj_coeffs = calJMT(current_d, target_d, duration);
                for (int i = 0; i < N_SAMPLES; i++) {
                        
                        double t = i * duration/N_SAMPLES;
                        double s_val = 0, d_val = 0;
                        for (int j = 0; j < s_traj_coeffs.size(); j++) {
                            s_val += this->s_traj_coeffs[j] * pow(t, j);
                            d_val += this->d_traj_coeffs[j] * pow(t, j);
                        }
                        s_traj.push_back(s_val);
                        d_traj.push_back(d_val);
                    }
                vector<vector<double>> possible_traj = {s_traj, d_traj};
                // double current_cost = calculate_total_cost(s_traj, d_traj, otherCarPredictions);
                double current_cost = calCost(s_traj, d_traj, otherCarPredictions);
                if (current_cost < best_cost) {
                        best_cost = current_cost;
                        best_frenet_traj = possible_traj;
                        best_traj_state = state;
                        best_target = target;
                    }
            }
        void calTrajWaypoints(vector<double> previous_path_x, vector<double> previous_path_y){
                next_x_vals.clear();
                next_y_vals.clear();
                coarse_s_traj.clear();
                coarse_x_traj.clear();
                coarse_y_traj.clear();
                interpolated_s_traj.clear();
                interpolated_x_traj.clear();
                interpolated_y_traj.clear();
                double prev_s = pos_s - s_dot * PATH_DT;                
                if (subpath_size >= 2) {
						coarse_s_traj.push_back(prev_s);
						coarse_x_traj.push_back(previous_path_x[subpath_size-2]);
						coarse_y_traj.push_back(previous_path_y[subpath_size-2]);
						coarse_s_traj.push_back(pos_s);
						coarse_x_traj.push_back(previous_path_x[subpath_size-1]);
						coarse_y_traj.push_back(previous_path_y[subpath_size-1]);
					} else {
						double prev_s = pos_s - 1;
						double prev_x = car_x - cos(deg2rad(car_yaw));
						double prev_y = car_y - sin(deg2rad(car_yaw));
						coarse_s_traj.push_back(prev_s);
						coarse_x_traj.push_back(prev_x);
						coarse_y_traj.push_back(prev_y);
						coarse_s_traj.push_back(pos_s);
						coarse_x_traj.push_back(car_x);
						coarse_y_traj.push_back(car_y);
					}
                // last two points of coarse trajectory, use target_d and current s + 30,60
                double target_s1 = pos_s + 30;
                double target_d1 = best_target[1][0];
                vector<double> target_xy1 = getXY(target_s1, target_d1, interpolated_waypoints_s, interpolated_waypoints_x, interpolated_waypoints_y);
                double target_x1 = target_xy1[0];
                double target_y1 = target_xy1[1];
                coarse_s_traj.push_back(target_s1);
                coarse_x_traj.push_back(target_x1);
                coarse_y_traj.push_back(target_y1);
                double target_s2 = target_s1 + 30;
                double target_d2 = target_d1;
                vector<double> target_xy2 = getXY(target_s2, target_d2, interpolated_waypoints_s, interpolated_waypoints_x, interpolated_waypoints_y);
                double target_x2 = target_xy2[0];
                double target_y2 = target_xy2[1];
                coarse_s_traj.push_back(target_s2);
                coarse_x_traj.push_back(target_x2);
                coarse_y_traj.push_back(target_y2);
                double target_s_dot = best_target[0][1];

                double current_s = pos_s;
                double current_v = s_dot;
                double current_a = s_ddot;
                cout << "current_v:" << current_v << endl;
                for (int i = 0; i < (NUM_PATH_POINTS - subpath_size); i++) {
                    if(current_v<SPEED_LIMIT && !too_close)
                        current_v += VELOCITY_INCREMENT_LIMIT;
                    else
                        current_v -= 0.8*VELOCITY_INCREMENT_LIMIT;
                        
                    // else{
                    //     current_v -= VELOCITY_INCREMENT_LIMIT;
                    //     if(current_v < follow_speed && !too_close)
                    //         current_v += VELOCITY_INCREMENT_LIMIT;
                    // }
                    current_s += current_v * PATH_DT;
                    interpolated_s_traj.push_back(current_s);
                }
                interpolated_x_traj = interpolate_points(coarse_s_traj, coarse_x_traj, interpolated_s_traj);
                interpolated_y_traj = interpolate_points(coarse_s_traj, coarse_y_traj, interpolated_s_traj);
                // cout << " is:" << interpolated_s_traj.size()
                //     << " ix:" << interpolated_x_traj.size()
                //     << " iy:" << interpolated_y_traj.size()
                //     << endl;
                for(int i = 0; i < subpath_size; i++) {
						next_x_vals.push_back(previous_path_x[i]);
						next_y_vals.push_back(previous_path_y[i]);
					} 
					// add xy points from newly generated path
					for (int i = 0; i < interpolated_x_traj.size(); i++) {
						//if (subpath_size == 0 && i == 0) continue; // maybe skip start position as a path point?
						next_x_vals.push_back(interpolated_x_traj[i]);
						next_y_vals.push_back(interpolated_y_traj[i]);
					} 
            }

        
        vector<double> calJMT(vector<double> start, vector<double> end, double T){
                MatrixXd a(3,3);
                double T2 =  T*T, 
                    T3 = T2*T, 
                    T4 = T3*T,
                    T5 = T4*T;
                a <<  T3,    T4,    T5, 
                    3*T2,  4*T3,  5*T4, 
                    6*T, 12*T2, 20*T3;
                MatrixXd aInv = a.inverse();
                
                VectorXd b(3);
                b << end[0] - (start[0] + start[1]*T + 0.5*start[2]*T2),
                    end[1] - (           start[1]   +     start[2]*T),
                    end[2] - (                            start[2]);
                VectorXd alpha = aInv * b;
                
                vector<double> output = {start[0], start[1], 0.5*start[2], alpha[0], alpha[1], alpha[2]};
                return output;
            }
        vector<double> getXvals(){return next_x_vals;}
        vector<double> getYvals(){return next_y_vals;}
        double calCost(vector<double> s_traj, vector<double> d_traj, vector<vector<vector<double>>> predictions){
                double total_cost = 0;
                double col = collision_cost(s_traj, d_traj, predictions) * COLLISION_COST_WEIGHT;
                double buf = buffer_cost(s_traj, d_traj, predictions) * BUFFER_COST_WEIGHT;
                double eff = efficiency_cost(s_traj) * EFFICIENCY_COST_WEIGHT;
                double ilb = in_lane_buffer_cost(s_traj, d_traj, predictions) * IN_LANE_BUFFER_COST_WEIGHT;
                double mjs = max_jerk_cost(s_traj) * MAX_JERK_COST_WEIGHT;
                double mjd = max_jerk_cost(d_traj) * MAX_JERK_COST_WEIGHT;
                double nml = not_middle_lane_cost(d_traj) * NOT_MIDDLE_LANE_COST_WEIGHT;
                // total_cost = buf + eff + mjs + mjd + nml;
                // total_cost = buf + eff + ilb + mjs + mjd + nml;
                total_cost = col + buf + eff + ilb + mjs + mjd + nml;
                // cout << " col:" << col << " buf:" << buf << " eff:" << eff << " ilb:" << ilb << endl;
                return total_cost;
            }
    private:
        double collision_cost(vector<double> s_traj, vector<double> d_traj, vector<vector<vector<double>>> predictions) {
                // Binary cost function which penalizes collisions.
                // double nearest = nearest_approach_to_any_vehicle_in_lane(s_traj, d_traj, predictions);
                double nearest = nearest_approach_to_any_vehicle(s_traj, d_traj, predictions);
                // cout << "nearest:" << nearest << endl;
                if (nearest < 2 * VEHICLE_RADIUS) {
                    return 1;
                } else {
                    return 0;
                }
            }
        double buffer_cost(vector<double> s_traj, vector<double> d_traj, vector<vector<vector<double>>> predictions) {
                // Penalizes getting close to other vehicles.
                double nearest = nearest_approach_to_any_vehicle(s_traj, d_traj, predictions);
                return logistic(2 * VEHICLE_RADIUS / nearest);
            }
        double nearest_approach_to_any_vehicle(vector<double> s_traj, vector<double> d_traj, vector<vector<vector<double>>> predictions){
                double closest = 999999;
                for (auto prediction : predictions) {
                    double current_dist = nearest_approach(s_traj, d_traj, prediction);
                    if (current_dist < closest) {
                        closest = current_dist;
                    }
                }
                return closest;
            }
        double nearest_approach(vector<double> s_traj, vector<double> d_traj, vector<vector<double>> prediction) {
                double closest = 999999;
                for (int i = 0; i < N_SAMPLES; i++) {
                    // double p_s = prediction[i][0];
                    // double p_d = prediction[i][1];
                    // double traj_s = s_traj[i];
                    // double traj_d = d_traj[i];
                    // vector<double> p_xy = getXY(p_s, p_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    // vector<double> traj_xy = getXY(traj_s, traj_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    // double current_dist = sqrt(pow(traj_xy[0] - p_xy[0], 2) + pow(traj_xy[1] - p_xy[1], 2));
                    double current_dist = sqrt(pow(s_traj[i] - prediction[i][0], 2) + pow(d_traj[i] - prediction[i][1], 2));                    
                    if (current_dist < closest) {
                        closest = current_dist;
                        // cout << "dist:" << closest << endl;
                    }
                }
                return closest;
            }
        double efficiency_cost(vector<double> s_traj) {
                // Rewards high average speeds.
                vector<double> s_dot_traj = velocities_for_trajectory(s_traj);
                double final_s_dot = s_dot_traj[s_dot_traj.size() - 1];
                return logistic((SPEED_LIMIT - final_s_dot) / SPEED_LIMIT);
            } 
        vector<double> velocities_for_trajectory(vector<double> traj) {
                vector<double> velocities;
                for (int i = 1; i < traj.size(); i++) {
                    velocities.push_back((traj[i] - traj[i-1]) / DT);
                }
                return velocities;
            }
        double in_lane_buffer_cost(vector<double> s_traj, vector<double> d_traj, vector<vector<vector<double>>> predictions) {
                // Penalizes getting close to other vehicles.
                double nearest = nearest_approach_to_any_vehicle_in_lane(s_traj, d_traj, predictions);
                // return logistic(2 * VEHICLE_RADIUS / nearest);
                if (nearest < 2 * VEHICLE_RADIUS) {
                    return 1;
                } else {
                    return 0;
                }
            }
        double nearest_approach_to_any_vehicle_in_lane(vector<double> s_traj, vector<double> d_traj, vector<vector<vector<double>>> predictions) {
                // Determines the nearest the vehicle comes to any other vehicle throughout a trajectory
                double closest = 999999;
                for (auto prediction : predictions) {
                double my_final_d = d_traj[d_traj.size() - 1];
                int my_lane = my_final_d / 4;
                vector<vector<double>> pred_traj = prediction;
                double pred_final_d = pred_traj[pred_traj.size() - 1][1];
                int pred_lane = pred_final_d / 4;
                if (my_lane == pred_lane) {
                        double current_dist = nearest_approach(s_traj, d_traj, prediction);
                        if (current_dist < closest && current_dist < 120) {
                            closest = current_dist;
                        }
                    }
                }
                return closest;
            }
        double max_jerk_cost(vector<double> s_traj) {
                // Penalize exceeding MAX_INSTANTANEOUS_JERK
                vector<double> s_dot_traj = velocities_for_trajectory(s_traj);
                vector<double> s_ddot_traj = velocities_for_trajectory(s_dot_traj);
                vector<double> s_dddot_traj = velocities_for_trajectory(s_ddot_traj);
                for (double s_dddot : s_dddot_traj) {
                    if (s_dddot > MAX_INSTANTANEOUS_JERK) {
                        return 1;
                    }
                }
                return 0;
            }
        double not_middle_lane_cost(vector<double> d_traj) {
                // penalize not shooting for middle lane (d = 6)
                double end_d = d_traj[d_traj.size()-1];
                return logistic(pow(end_d-6, 2));
            }
    
    
    private:
            vector<double> next_x_vals;
          	vector<double> next_y_vals;
            vector<double> coarse_waypoints_s;
            vector<double> coarse_waypoints_x;
            vector<double> coarse_waypoints_y;
            vector<double> coarse_waypoints_dx;
            vector<double> coarse_waypoints_dy;

            vector<double> interpolated_waypoints_s;
            vector<double> interpolated_waypoints_x;
            vector<double> interpolated_waypoints_y;
            vector<double> interpolated_waypoints_dx;
            vector<double> interpolated_waypoints_dy;

            vector<double> map_waypoints_s;
            vector<double> map_waypoints_x;
            vector<double> map_waypoints_y;
            vector<double> map_waypoints_dx;
            vector<double> map_waypoints_dy;

            double car_x, car_y, car_yaw, car_s, car_d, car_speed;
            double pos_s, s_dot, s_ddot;
            double pos_d, d_dot, d_ddot;
            double traj_start_time; 
            double duration;
            double best_cost;   
            double follow_speed;
            string best_traj_state;
                         
            int target_lane, current_lane;
            int subpath_size;
            int num_waypoints;            

            bool car_to_left, car_to_right, car_just_ahead;
            bool too_close;

            vector<vector<vector<double>>> otherCarPredictions;
            vector<vector<double>> best_frenet_traj, best_target;
            vector<vector<double>> sensor_fusion;
            vector<double> s_traj_coeffs;
            vector<double> d_traj_coeffs;

            vector<double> coarse_s_traj;
            vector<double> coarse_x_traj;
            vector<double> coarse_y_traj;
            vector<double> interpolated_s_traj; 
			vector<double> interpolated_x_traj;
            vector<double> interpolated_y_traj;

            vector<vector<double>> carsInLane0;
            vector<vector<double>> carsInLane1;
            vector<vector<double>> carsInLane2;
            
            
            vector<string> available_states;
            
};
#endif
#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y) {

  double closestLen = 100000; //large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x, y, map_x, map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }

  }

  return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y) {

  int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y - y), (map_x - x));

  double angle = fabs(theta - heading);
  angle = min(2 * pi() - angle, angle);

  if (angle > pi() / 4) {
    closestWaypoint++;
    if (closestWaypoint == maps_x.size()) {
      closestWaypoint = 0;
    }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y) {
  int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

  int prev_wp;
  prev_wp = next_wp - 1;
  if (next_wp == 0) {
    prev_wp = maps_x.size() - 1;
  }

  double n_x = maps_x[next_wp] - maps_x[prev_wp];
  double n_y = maps_y[next_wp] - maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
  double proj_x = proj_norm * n_x;
  double proj_y = proj_norm * n_y;

  double frenet_d = distance(x_x, x_y, proj_x, proj_y);

  //see if d value is positive or negative by comparing it to a center point

  double center_x = 1000 - maps_x[prev_wp];
  double center_y = 2000 - maps_y[prev_wp];
  double centerToPos = distance(center_x, center_y, x_x, x_y);
  double centerToRef = distance(center_x, center_y, proj_x, proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
  }

  frenet_s += distance(0, 0, proj_x, proj_y);

  return {frenet_s, frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s,
                     double d,
                     const vector<double> &maps_s,
                     const vector<double> &maps_x,
                     const vector<double> &maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
    prev_wp++;
  }

  int wp2 = (prev_wp + 1) % maps_x.size();

  double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s - maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
  double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

  double perp_heading = heading - pi() / 2;

  double x = seg_x + d * cos(perp_heading);
  double y = seg_y + d * sin(perp_heading);

  return {x, y};

}

double getLaneFrenetD(int lane_id) {
  return 2 + 4 * lane_id;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy](uWS::WebSocket<
      uWS::SERVER> ws,
                                                                                                           char *data,
                                                                                                           size_t length,
                                                                                                           uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];
          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          json msgJson;

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          const double speed_limit = 49;

          // Do decision planning first, then do path generation
          double eventual_target_speed;

          // Use perception result to do decision planning

          // Set object lane
          int target_lane_id = 1;

          // Set target speed
          double target_speed;
          double speed_step = 2;
          if (car_speed + speed_step < eventual_target_speed) {
            target_speed = car_speed + speed_step;
          } else {
            target_speed = car_speed - speed_step;
          }

          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

          // Temporary points to work with for calculating the fitting curve
          vector<double> ptsx;  // In world coodinates
          vector<double> ptsy;  // In world coordinates

          // Now we are preparing the points that are the input for the curve fitting
          // We plan to take the last 2 points from the previous planning points, and then add several target points
          // to create a list of points to fit

          int prev_path_size = previous_path_x.size();
          car_yaw = deg2rad(car_yaw);

          if (prev_path_size < 2) {
            // If the last planning path has less than 2 points, we will just make up the these 2 points now
            // We make up 2 points: The car's current location, and a imaginary previous location as if the car was
            // traveling in straight line in the last step
            double imaginary_last_step_duration = 1.0;
            double current_car_x_w = car_x;
            double current_car_y_w = car_y;
            double last_car_x_w = car_x - cos(car_yaw) * imaginary_last_step_duration;
            double last_car_y_w = car_y - sin(car_yaw) * imaginary_last_step_duration;
            ptsx.push_back(last_car_x_w);
            ptsx.push_back(current_car_x_w);
            ptsy.push_back(last_car_y_w);
            ptsy.push_back(current_car_y_w);
          } else {
            // Otherwise, we will just pick up the last two points from previous path planning. We do this because we
            // want to ensure the new planned path can smoothly connect with the previous planned path
            ptsx.push_back(previous_path_x[prev_path_size - 2]);
            ptsx.push_back(previous_path_x[prev_path_size - 1]);
            ptsy.push_back(previous_path_y[prev_path_size - 2]);
            ptsy.push_back(previous_path_y[prev_path_size - 1]);
          }

          // Also, the first 2 points of the created list is going to be used to determine the local coordination system.
          // The yaw of the vector (ref_yaw) created by the first 2 points decide the rotation from world coordinates,
          // and the location of the first point (ref_x and ref_y) decides the translation from world coordinates
          double ref_x = ptsx[1];
          double ref_y = ptsy[1];
          double ref_yaw = atan2(ptsy[1] - ptsy[0], ptsx[1] - ptsx[0]);
          auto ref_sd = getFrenet(ref_x, ref_y, ref_yaw, map_waypoints_x, map_waypoints_y);
          double ref_s = ref_sd[0];
          double ref_d = ref_sd[1];

          // Add more target points into the list for fitting
          vector<double> next_waypoint_1 =
              getXY(car_s + 30, getLaneFrenetD(target_lane_id), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_waypoint_2 =
              getXY(car_s + 60, getLaneFrenetD(target_lane_id), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_waypoint_3 =
              getXY(car_s + 90, getLaneFrenetD(target_lane_id), map_waypoints_s, map_waypoints_x, map_waypoints_y);

          ptsx.push_back(next_waypoint_1[0]);
          ptsx.push_back(next_waypoint_2[0]);
          ptsx.push_back(next_waypoint_3[0]);

          ptsy.push_back(next_waypoint_1[1]);
          ptsy.push_back(next_waypoint_2[1]);
          ptsy.push_back(next_waypoint_3[1]);

//          std::cout << "prev_count " << previous_path_x.size() << " " << previous_path_y.size() << std::endl;

          // Translate the points from world coordinates to local coordinates (per ref_x, ref_y, ref_yaw) so we can fit them more easily
          vector<double> ptsx_l;
          vector<double> ptsy_l;
          for (int i = 0; i < ptsx.size(); ++i) {
            double x_l = (ptsx[i] - ref_x) * cos(ref_yaw) + (ptsy[i] - ref_y) * sin(ref_yaw);
            double y_l = -(ptsx[i] - ref_x) * sin(ref_yaw) + (ptsy[i] - ref_y) * cos(ref_yaw);

            ptsx_l.push_back(x_l);
            ptsy_l.push_back(y_l);
          }

          // Fit with a spline
          tk::spline s;
          s.set_points(ptsx_l, ptsy_l);

          // Generate output points
          // First, add back remaining points from last planning
          for (int i = 0; i < previous_path_x.size(); ++i) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          // Then sample some of the points on the fitted spline as new waypoints for the new planning path
          // We design a way such that, every time we sample a point, it's one step further away from ref_x and ref_y
          double target_x_l = 30;
          double target_y_l = s(target_x_l);
          double target_distance = sqrt(target_x_l * target_x_l + target_y_l * target_y_l);

          double step_duration = 0.02;
          double step_length = target_speed * 0.447 * step_duration;

          // Assume we are heading straight to the target_x_l and tart_y_l, how much steps we need?
          // of course this is just an estimate because we probably won't go straight there
          int step_count = static_cast<int>(target_distance / step_length);

          // Now we know for each sampling, how much incremental local x we want
          double step_x = target_x_l / step_count;

          // Start the sampling and add sampled points
          const int new_points_count = 50 - previous_path_x.size();
          for (int i = 0; i < new_points_count; ++i) {
            double sample_x_l = (i + 1) * step_x;
            double sample_y_l = s(sample_x_l);

            // Translate from local coordinates to world coordinates
            double x_w = ref_x + cos(ref_yaw) * sample_x_l - sin(ref_yaw) * sample_y_l;
            double y_w = ref_y + sin(ref_yaw) * sample_x_l + cos(ref_yaw) * sample_y_l;

            next_x_vals.push_back(x_w);
            next_y_vals.push_back(y_w);
          }

//          std::cout << "count " << next_x_vals.size() << " " << next_y_vals.size() << std::endl;

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\"," + msgJson.dump() + "]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

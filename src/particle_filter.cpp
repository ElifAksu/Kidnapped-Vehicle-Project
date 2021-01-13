/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {

  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  std::default_random_engine gen;
  num_particles = 100;  // TODO: Set the number of particles

  double std_x= std[0];
  double std_y= std[1];
  double std_theta = std[2];

  std::normal_distribution<double> gaus_x(x, std_x);
  std::normal_distribution<double> gaus_y(y, std_y);
  std::normal_distribution<double> gaus_theta(theta, std_theta);


  particles = vector<Particle>(num_particles);
  Particle particle;
  for (int i=0; i<num_particles; ++i){
      particle.x = gaus_x(gen);
      particle.y = gaus_y(gen);
      particle.theta = gaus_theta(gen);
      particle.weight = 1.0; // equal weights !
      particles[i] = particle;
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {

  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  std::default_random_engine gen;
  double std_x= std_pos[0];
  double std_y= std_pos[1];
  double std_theta= std_pos[2];


  double x_pred;
  double y_pred;
  double theta_pred;

  for (int i=0; i<num_particles; ++i){
      if (fabs(yaw_rate)<0.0001){ // yaw rate can be zero
          x_pred = particles[i].x + velocity * delta_t * cos(particles[i].theta);
          y_pred = particles[i].y + velocity * delta_t * sin(particles[i].theta);
      }
      else { 
          x_pred = particles[i].x + (velocity/yaw_rate) * ( sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta) );
          y_pred = particles[i].y + (velocity/yaw_rate) * ( -cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta) );
      }
      theta_pred = particles[i].theta + yaw_rate * delta_t;

     
      std::normal_distribution<double> gaus_x(x_pred, std_x);
      std::normal_distribution<double> gaus_y(y_pred, std_y);
      std::normal_distribution<double> gaus_theta(theta_pred, std_theta);

      particles[i].x = gaus_x(gen);
      particles[i].y = gaus_y(gen);
      particles[i].theta = gaus_theta(gen);
      }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  for(auto& obs: observations){
    double min_dist_value = 1000.0; 

    for(const auto& pred: predicted){
      double d = dist(obs.x, obs.y, pred.x, pred.y); //from helper.
      if( min_dist_value > d){
        min_dist_value = d;
        obs.id = pred.id;
      }
    }
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a multi-variate Gaussian
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  for(auto& p: particles){
    p.weight = 1.0;

    // Prediction
    vector<LandmarkObs> predictions;
    for(const auto& landmark: map_landmarks.landmark_list){
      double distance = dist(p.x, p.y, landmark.x_f, landmark.y_f);
      if( distance < sensor_range){ // prediction- sensor range check
        predictions.push_back(LandmarkObs{landmark.id_i, landmark.x_f, landmark.y_f});
      }
    }

    // transform observationss
    vector<LandmarkObs> observations_transform;

    for(const auto& observ: observations){
      LandmarkObs tmp;
      tmp.x = observ.x * cos(p.theta) - observ.y * sin(p.theta) + p.x;
      tmp.y = observ.x * sin(p.theta) + observ.y * cos(p.theta) + p.y;
      observations_transform.push_back(tmp);
    }

    // Nearist Neighbour
    dataAssociation(predictions, observations_transform);

    // multi-variate Gaussian :
    for(const auto& obs_tr: observations_transform){

      Map::single_landmark_s landmark = map_landmarks.landmark_list.at(obs_tr.id-1);
      double x_e= pow(obs_tr.x - landmark.x_f, 2) / (2 * pow(std_landmark[0], 2)); //error terms
      double y_e = pow(obs_tr.y - landmark.y_f, 2) / (2 * pow(std_landmark[1], 2)); //error terms
      double normalizer =(2 * M_PI * std_landmark[0] * std_landmark[1]);
      double w = exp(-(x_e + y_e)) /normalizer ;
      p.weight *=  w;
    }

    weights.push_back(p.weight);
  }
}


void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  vector<Particle> new_particles;

  vector<double> weights;
  for(int i=0; i<num_particles; i++){
    weights.push_back(particles[i].weight);
  }
  std::uniform_real_distribution<double> unifrom_dist(0.0, 1.0);
  int index   = int(unifrom_dist(gen) * num_particles);
  double beta = 0.0;
  double mw   = *max_element(weights.begin(), weights.end());
  for(int i=0; i<num_particles; i++){
    beta += unifrom_dist(gen) * 2.0 * mw;
    while(beta > weights[index]){
      beta -= weights[index];
      index = (index+1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }
  particles = new_particles;
  weights.clear();
}


void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  //Clear the previous associations
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
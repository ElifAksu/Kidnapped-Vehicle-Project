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
  std::cout<< "first measurement x " << x << " y: " << y<<std::endl;
  num_particles = 10;  // TODO: Set the number of particles
  std::default_random_engine rd{};
  std::mt19937 gen{rd()};
  particles.resize(num_particles);
    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
  std::normal_distribution<> gaus_x{x,std[0]};
  std::normal_distribution<> gaus_y{y,std[1]};
  std::normal_distribution<> gaus_theta{theta,std[2]};
  for(int i=0; i<num_particles;i++)
  {
   particles[i].x=gaus_x(gen);
   particles[i].y=gaus_y(gen);
   particles[i].theta=gaus_theta(gen);
   particles[i].weight=(1/num_particles);
   std::cout<< "particle " <<i << " x value: " << particles[i].x<<std::endl;
  }

is_initialized=true;
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
  ;
  for(int i=0; i<num_particles;i++)
  {
    double pred_x;
    double pred_y;
    double pred_theta;
    std::default_random_engine gen;
    if (fabs(yaw_rate) > 0.0001)
    {
    
    double newtheta=particles[i].theta+delta_t*yaw_rate;
    pred_x=particles[i].x + (velocity/yaw_rate)*(sin(newtheta)-sin(particles[i].theta));
    pred_y=particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) -cos(newtheta));
    pred_theta=particles[i].theta+newtheta;
    }
    else
    {
    pred_x=particles[i].x + velocity*delta_t*cos(particles[i].theta);
    pred_y=particles[i].y + velocity*delta_t*sin(particles[i].theta);
    pred_theta=particles[i].theta;
    }
    std::normal_distribution<> gaus_x{pred_x,std_pos[0]};
    std::normal_distribution<> gaus_y{pred_y,std_pos[1]};
    std::normal_distribution<> gaus_theta{pred_theta,std_pos[2]};
    particles[i].x=gaus_x(gen);
    particles[i].y=gaus_y(gen);
    particles[i].theta=gaus_theta(gen);
   
   std::cout<< "predicted particle " <<i << " x value: " << particles[i].x<<std::endl;
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations,Particle& particle ) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  int nearist;
  vector<int> associations; 
  vector<double> sense_x; 
  vector<double> sense_y;
  for(int i ; i<predicted.size(); i++)
  {
    double pred_x=predicted[i].x;
    double pred_y=predicted[i].y;
    double max_dist=1000;
    for(int j=0; j<observations.size();j++)
    {
      double obser_x = observations[j].x;
      double obser_y = observations[j].y;
      double distance = dist(pred_x,pred_y, obser_x,obser_y);
      if(distance < max_dist)
      {
        nearist=j;
      }
    }
    associations.push_back(observations[nearist].id);
    sense_x.push_back(observations[nearist].x);
    sense_y.push_back(observations[nearist].y);


  }
  SetAssociations(particle, associations, sense_x, sense_y);
 associations.clear();
 sense_x.clear();
 sense_y.clear();

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
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
   * 
   */
    for(int i=0; i<num_particles;i++)
    {
      vector<LandmarkObs> predicted_observations; //Map frame -- sensor range, predicted position
      for(int j=0; j< map_landmarks.landmark_list.size(); j++)
      {
        double diff_x= fabs(particles[i].x-map_landmarks.landmark_list[j].x_f);
        double diff_y= fabs(particles[i].y-map_landmarks.landmark_list[j].y_f);
        if(diff_x < sensor_range && diff_y<sensor_range)
        {
          LandmarkObs obs{map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f};
          predicted_observations.push_back(obs);
        }
      }
      vector<LandmarkObs> transformed_observations; // Observation in Map frame 
      for(int j=0; j< observations.size(); j++)
      {
        double map_obs_x= particles[i].x + (observations[j].x *cos(particles[i].theta))-(observations[j].y *sin(particles[i].theta));
        double map_obs_y= particles[i].y + (observations[j].x *sin(particles[i].theta))+(observations[j].y *cos(particles[i].theta));

        LandmarkObs tr_obs{observations[j].id, map_obs_x,map_obs_y};
        transformed_observations.push_back(tr_obs);
      }

      dataAssociation(predicted_observations, transformed_observations, particles[i]);
    }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

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
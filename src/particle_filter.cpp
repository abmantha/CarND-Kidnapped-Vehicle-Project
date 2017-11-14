/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    num_particles = 10; // some empirical value

    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]); 

    particles.resize(num_particles);
    for (auto i = 0; i < particles.size(); ++i)
    {
        Particle p;
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = 1;
        particles.push_back(p);
    }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;


    for (auto i = 0; i < particles.size(); ++i)
    {
        Particle p = particles[i];

        normal_distribution<double> dist_x(p.x, std_pos[0]);
        normal_distribution<double> dist_y(p.y, std_pos[1]);
        normal_distribution<double> dist_yaw(p.theta, std_pos[2]); 

        p.x = (p.x + (velocity / yaw_rate) * (std::sin(p.theta + yaw_rate*delta_t) - std::sin(p.theta))) + dist_x(gen);
        p.y = (p.y + (velocity / yaw_rate) * (std::cos(p.theta) - std::cos(p.theta + yaw_rate*delta_t))) + dist_y(gen);
        p.theta = (p.theta + yaw_rate*delta_t) + dist_yaw(gen);
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    // predicted -- prediction measurements between one particular particle and all landmarks within range
    // observations -- actual landmark measurements gathered from lidar

    // assign each sensor observation the map landmark id associated with it

    // for each predicted value
    //      store the min value to be the current value
    //      for each observation value
    //          calculate the distance between the prediction and observation
    //          if the difference is greater than the min
    //              set the min to be the difference

    for (auto p = predicted.begin(); p != predicted.end(); ++p)
    {
        double min_distance = std::sqrt(std::pow(p->x_f, 2) + std::pow(p->y_f, 2));
        int map_landmark_id;
        for (auto o = observations.begin(); o != observations.end(); ++o)
        {
            double curr_distance = std::sqrt(std::pow(o->x_f - p->x_f, 2) + std::pow(o->y_f - p->y_f, 2));
            if (curr_distance < min_distance)
            {
                min_distance = curr_distance;
                map_landmark_id = o.id;
            }
        }
        p.id = map_landmark_id;
    }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    double std_x = std_landmark[0];
    double std_y = std_landmark[1];
    double gaussian_norm = 1 / (2 * M_PI * std_x * std_y);

    for (auto p = particles.begin(); p != particles.end(); ++p) 
    {
        double weight_product = 1.0;
        for (auto o = observations.begin(); o != observations.end(); ++o)
        {
            double transform_x = p.x + (std::cos(p.theta) * o.x_f) - (std::sin(p.theta) * o.y_f);
            double transform_y = p.y + (std::sin(p.theta) * o.x_f) + (std::cos(p.theta) * o.y_f);

            double min_distance = 1000000;
            int min_landmark_index;

            for (auto l = map_landmarks.landmark_list.begin(); l != map_landmarks.landmark_list.end(); ++l)
            {
                double diff = std::sqrt(std::pow(transform_x - l.x_f, 2) + std::pow(transform_y - l.y_f, 2));
                if (diff < min_distance)
                {
                    min_distance = diff;
                    min_landmark_index = l.id_i;
                    p.sense_x.push_back(transform_x);
                    p.sense_y.push_back(transform_y);
                }
            }

            p.associations.push_back(min_landmark_index);

            // calculate weight for particle
            double mu_x = map_landmarks.landmark_list[min_landmark_index].x_f;
            double mu_y = map_landmarks.landmark_list[min_landmark_index].y_f;
            double exponent = (std::pow(transform_x - mu_x, 2) / (2 * std::pow(std_x, 2)) + std::pow(transform_y - mu_y, 2) / (2 * std::pow(std_y, 2));
            weight_product *= gaussian_norm * std::exp(-exponent);
        }
        p.weight = weight_product;
        weights.at(p.id) = weight_product;
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    int index = int(std::rand() % particles.size());
    double beta = 0.0;
    double max_weight = *std::max_element(weights.begin(), weights.end());

    std::uniform_real_distribution<> dist(0.0, 1.0);

    std::vector<Particle> resampled_particles;
    for (auto i = 0; i < particles.size(); ++i)
    {
        double r = ;
        beta += dist(gen) * 2.0 * max_weight;
        while (beta > weights[index])
        {
            beta -= weights[index];
            index = (index + 1) % particles.size();
        }
        resampled_particles.push_back(particles[index]);
    }
    particles = resampled_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
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

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

#ifndef __MOVINGSHAPEGEN_HPP
#include<Eigen/Dense>
#include<dmeans/utils>
#include<vector>
#include "movingdatagen.hpp"

class MovingShapeDataGenerator : public MovingDataGenerator{
	typedef Eigen::Vector2d V2d;
	enum Shape{
		LINE,
		CIRCLE,
		SQUARE,
		TRIANGLE
	};
	protected:
	std::vector<V2d> centers;
	std::vector<Shape> shapes;
	std::vector<double> rots;
	double motionStdDev, clusterStdDev, radius;

	public:
		MovingShapeDataGenerator(dmeans::Config cfg) : MovingDataGenerator(cfg){
			radius = cfg.get("radius", dmeans::Config::OPTIONAL, 0.1);
			motionStdDev = cfg.get("motionStdDev", dmeans::Config::OPTIONAL, 0.05);
			clusterStdDev = cfg.get("clusterStdDev", dmeans::Config::OPTIONAL, 0.05);

			std::uniform_real_distribution<double> uniformDist01(0, 1);
			for (uint32_t i = 0; i < this->alive.size(); i++){
					V2d newCenter;
					newCenter(0) = uniformDist01(dmeans::RNG::get());
					newCenter(1) = uniformDist01(dmeans::RNG::get());
					centers.push_back(newCenter);
					rots.push_back(2*M_PI*uniformDist01(dmeans::RNG::get()));
					shapes.push_back((int)(4.0*uniformDist01(dmeans::RNG::get())));
			}
		}

		virtual ~MovingShapeDataGenerator(){
		}

		void transitionCluster(int k){
			std::uniform_real_distribution<double> uniformDistAng(0, 2*M_PI);
			std::normal_distribution<double> transitionDistRadial(0, motionStdDev);

			double steplen = transitionDistRadial(dmeans::RNG::get());
			//angle sampled from uniform
			double stepang = uniformDistAng(dmeans::RNG::get());
			centers[k](0) += steplen*cos(stepang);
			centers[k](1) += steplen*sin(stepang);

			double rotang = transitionDistRadial(dmeans::RNG::get());
			rots[k] += rotang;

			cout << "Cluster " << k << " moved to " << centers[k].transpose() << endl;
		}

		void createCluster(){
				std::uniform_real_distribution<double> uniformDist01(0, 1);
				V2d newCenter;
				newCenter(0) = uniformDist01(dmeans::RNG::get());
				newCenter(1) = uniformDist01(dmeans::RNG::get());
				centers.push_back(newCenter);
				shapes.push_back((int)(4.0*uniformDist01(dmeans::RNG::get())));
				rots.push_back(2*M_PI*uniformDist01(dmeans::RNG::get()));
				cout << "Cluster " << centers.size()-1 << " was created at " << centers.back().transpose() << endl;
		}

		std::vector<V2d> genData(int j, int n) const{
		    uniform_real_distribution<double> uniformDist01(0, 1);
			uniform_real_distribution<double> uniformDistAng(0, 2*M_PI);
			normal_distribution<double> likelihoodDistRadial(0, clusterStdDev);
			std::vector<V2d> out;
			switch(shapes[j]){
				case SQUARE:
					for (int k = 0; k < n; k++){
						V2d noise;
						double len = likelihoodDistRadial(dmeans::RNG::get());
						double ang = uniformDistAng(dmeans::RNG::get());
						noise(0) += len*cos(ang);
						noise(1) += len*sin(ang);
						double  ang2 = uniformDistAng(dmeans::RNG::get());
						V2d newData = centers[j] + noise;
						if (ang2 < 45.0*M_PI/180.0){
							newData(0) += radius;
							newData(1) += radius*tan(rots[j]);
						} else if (ang2 < 135.0*M_PI/180.0){
							ang2 -= M_PI/2.0;
							newData(0) -= radius*tan(rots[j]);
							newData(1) += radius;
						} else if (ang2 < 225.0*M_PI/180.0){
							ang2 -= M_PI;
							newData(0) -= radius;
							newData(1) -= radius*tan(rots[j]);
						} else if (ang2 < 315.0*M_PI/180.0){
							ang2 -= 3.0*M_PI/2.0;
							newData(0) += radius*tan(rots[j]);
							newData(1) -= radius;
						} else {
							newData(0) += radius;
							newData(1) += radius*tan(rots[j]);
						}
						out.push_back(newData);
					}
					break;
				case LINE:
					for (int k = 0; k < n; k++){
						V2d noise;
						double len = likelihoodDistRadial(dmeans::RNG::get());
						double ang = uniformDistAng(dmeans::RNG::get());
						noise(0) += len*cos(ang);
						noise(1) += len*sin(ang);
						V2d newData;
						newData = centers[j] + noise;
						newData(0) += radius*(2.0*uniformDist01(dmeans::RNG::get()) - 1.0)
						out.push_back(newData);
					}
					break;
				case TRIANGLE:
					for (int k = 0; k < n; k++){
						V2d noise, newData;
						double len = likelihoodDistRadial(dmeans::RNG::get());
						double ang = uniformDistAng(dmeans::RNG::get());
						noise(0) += len*cos(ang);
						noise(1) += len*sin(ang);
						double sd = uniformDist01(dmeans::RNG::get());
						newData = centers[j];
						newData(1) -= radius/sqrt(3.0);
						if(sd < .333){
							newData(0) += radius*(2.0*uniformDist01(dmeans::RNG::get()) - 1.0);
						} else if (sd < .666){
							double t = radius*uniformDist01(dmeans::RNG::get());
							newData(0) += -radius + t;
							newData(1) += t*tan(M_PI/3.0);
						} else {
							double t = radius*uniformDist01(dmeans::RNG::get());
							newData(0) += radius -t;
							newData(1) += t*tan(M_PI/3.0);
						}
						out.push_back(newData);
					}
					break;
				case CIRCLE:
					for (int k = 0; k < n; k++){
						V2d noise, newData;
						double len = likelihoodDistRadial(dmeans::RNG::get());
						double ang = uniformDistAng(dmeans::RNG::get());
						noise(0) += len*cos(ang);
						noise(1) += len*sin(ang);
						out.push_back(newData);

						V2d newData = centers[j];
						double len = radius+likelihoodDistRadial(dmeans::RNG::get());
						double ang = uniformDistAng(dmeans::RNG::get());
						newData(0) += len*cos(ang);
						newData(1) += len*sin(ang);
						out.push_back(newData);
					}
					break;
				default:
					cout << "Error in MovingShapeDataGenerator: Shape unknown" << endl;
					exit(0);
					break;
			}
			//rotate the data
			for(int k = 0; k < n; k++){
				V2d r;
				r(0) = out[k](0)*cos(rots[j])-out[k](1)*sin(rots[j]);
				r(1) = out[k](0)*sin(rots[j])+out[k](1)*cos(rots[j]);
				out[k] = r;
			}
			return out;
		}
};
#define __MOVINGSHAPEGEN_HPP
#endif /* __MOVINGSHAPEGEN_HPP */

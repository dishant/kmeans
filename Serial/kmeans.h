/*
 * kmeans.h
 *
 *  Created on: Nov 5, 2012
 *      Author: dishant
 */

#ifndef KMEANS_H_
#define KMEANS_H_

/*
 * k-means.cpp
 *
 *  Created on: Nov 5, 2012
 *      Author: dishant
 */


#include <iostream>
#include <vector>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <math.h>

using namespace std;


class Point{
      vector<double> coordinates;
      int clusterID;
      double cluster_distance;

      public:

      Point(vector<double>& coord)
      {
             coordinates = coord;
             clusterID = -1;
             cluster_distance = std::numeric_limits<double>::max();
      }

      void add(Point p1)
      {
           vector<double> p1_coordinates = p1.get_coordinates();
           for(int i=0;i<coordinates.size();i++)
           {
                   coordinates[i] += p1_coordinates[i];
           }
      }
      void set_cluster_id(int id)
      {
             clusterID=id;
      }
      void divide_by(int p)
            {
                   for(int i=0;i<coordinates.size();i++)
                   {
                           coordinates[i]/=p;
                   }
            }

      vector<double> get_coordinates()
      {
             return coordinates;
      }

      int get_dimensions()
      {
          return coordinates.size();
       }
      int get_cluster_id()
      {
             return clusterID;
      }

      void set_cluster_distance(double min_dist)
      {
    	  cluster_distance =  min_dist;
      }

      double get_cluster_distance()
		{
		  return cluster_distance;
		}

};

class Cluster{
      vector<double> coordinates;
      //vector<Point> points_in_cluster;
      // No. of points in the cluster

      public:

      Cluster(vector<double> coord)
      {
             coordinates=coord;
      }

      double calculate_dist(Point p)
      {
             double dist=0;
             vector<double> p_coordinates=p.get_coordinates();

             for(int i=0;i<coordinates.size();i++)
             {
                     dist += (coordinates[i] - p_coordinates[i])*(coordinates[i] - p_coordinates[i]);
             }
             return sqrt(dist);
      }


      void set_coordinates(vector<double> coord)
      {
             coordinates =coord;
      }

      vector<double> get_coordinates()
      {
             return coordinates ;
      }

};


bool greater_than(Point p1, Point p2)
    {
           vector<double> p1_coordinates = p1.get_coordinates();
           vector<double> coordinates = p2.get_coordinates();

           for(int i=0; i<coordinates.size();i++)
           {
                   if(coordinates[i] > p1_coordinates[i])
                   {
                                     return true;
                                     }
                   else if(coordinates[i] < p1_coordinates[i])
                    { return false;
                    }
           }
           return false;
    }

void set_min_cluster(vector<Cluster> clusters, Point& p)
{
    double min_dist=std::numeric_limits<double>::max();
    int index=-1;
    int n = clusters.size();

    for(int i=0;i<n;i++)
    {
        double dist=clusters[i].calculate_dist(p);
       if( dist < min_dist)
       {
           index=i;
           min_dist=dist;
       }
    }
    p.set_cluster_distance(min_dist);
    p.set_cluster_id(index);

}

int compute_closest_cluster(vector<Cluster> input_cluster,Cluster output_cluster, double& min_distance)
{
	int k = input_cluster.size();
	int closest_cluster = -1;

	for(int j=0;j<k;j++)
	{
		vector<double> s = input_cluster[j].get_coordinates();
		Point* p = new Point(s);
		double temp = output_cluster.calculate_dist(*p);
		if(min_distance > temp)
		{
			min_distance = temp;
			closest_cluster = j;
		}
	}

	return closest_cluster;

}

double check_accuracy(vector<Cluster> input_cluster, vector<Cluster> output_cluster,double x)
{
	int k = input_cluster.size();
	int* closest_cluster = new int[k];
	double* min_distance = new double[k];
	double allowed_distance = x;
	int correct_classifications = 0;


	for(int i=0;i<k;i++)
	{
		min_distance[i] = std::numeric_limits<double>::max();
		closest_cluster[i] = compute_closest_cluster(input_cluster,output_cluster[i],min_distance[i]);
		if(min_distance[i] < allowed_distance)
		{
			correct_classifications++;
		}
	}

return 100*correct_classifications/k;
}

vector<Cluster> kmeans(int k, vector<Point> points)
{

				sort(points.begin(),points.end(),greater_than);

                int n=points.size();
                int step=n/k;
                int dimension = points[0].get_dimensions();

                vector<Cluster> clusters;

                //Create the clusters
                for(int j=0; j<step*k;j=j+step)
                {

                        vector<double> sum(dimension,0);
                        for(int r=j; r<j+step;r++)
                        {

                                for(int z=0;z<dimension;z++)
                                {
                                        sum[z] += (points[r].get_coordinates())[z];

                                }


                        }
                        for(int z=0;z<dimension;z++)
						{
								sum[z] = sum[z]/step;
						}

                        clusters.push_back(Cluster(sum));
                }
                int t=0;
                while(t++<5)
                {
                	/*
                	cout<<"Printing current clusters\n";
                	for(int i=0; i<clusters.size(); i++)
					{
						vector<double> d = clusters[i].get_coordinates();
						cout<<"Cluster"<<i<<" ("<<d[0]<<","<<d[1]<<')'<<endl;
					}
					*/

                    //Associate the points to the clusters
                    for(int j=0;j<n;j++)
                    {
                          set_min_cluster(clusters,points[j]);

                    }


                    //Move the clusters according to the clusters and associated points
                    vector<double> temp(dimension,0);
                    vector<Point> dist_p_cluster(k,temp);
                    vector<int> no_of_points_per_cluster(k,0);

                    for(int j=0;j<n;j++)
                    {
                         (dist_p_cluster[points[j].get_cluster_id()]).add(points[j]);
                         no_of_points_per_cluster[points[j].get_cluster_id()]++;
                    }

                    for(int j=0;j<k;j++)
                    {
                    	if(no_of_points_per_cluster[j] > 0)
                    	{
							dist_p_cluster[j].divide_by(no_of_points_per_cluster[j]);
							clusters[j].set_coordinates(dist_p_cluster[j].get_coordinates());
                    	}
                    }

                    double max_distance = 0;
                    int farthest_point_index = -1;
                    for(int j=0; j<n; j++)
                    {
                    	if(points[j].get_cluster_distance() > max_distance)
                    	{
                    		max_distance = points[j].get_cluster_distance();
                    		farthest_point_index = j;
                    	}
                    }

                    for(int j=0;j<k;j++)
                    {
                    	if(no_of_points_per_cluster[j] == 0)
                    	{
                    		clusters[j].set_coordinates(points[farthest_point_index].get_coordinates());
                    		break;
                    	}
                    }
                }

                return clusters;
}


#endif /* KMEANS_H_ */


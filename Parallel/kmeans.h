#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <math.h>

/* Range of input space per dimension*/
#define MIN_VALUE 0
#define MAX_VALUE 200000

using namespace std;
class Point{
      vector<double> coordinates;
      int clusterID;
      double cluster_distance;

      public:

      Point(vector<double>& coord)
      {
             coordinates = coord;
             clusterID = 0;
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
      void set_coordinates(vector<double> coord)
      {
		   coordinates =coord;
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

void set_min_cluster(vector<Point> clusters, Point& p)
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

int compute_closest_cluster(vector<Point> input_cluster,Point output_cluster, double& min_distance)
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

double check_accuracy(vector<Point> input_cluster, vector<Point> output_cluster,double x)
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


/**
 * @arg1 dimensions of the input dataset
 * @arg2 clusters needed
 * @arg3 total number of points in the dataset > size of dataset
 */
vector<Point> get_dataset(int dimensions, int k, int total_points, vector<Point>& clusters, double& x)
{
	//vector<Point> clusters;
	vector<Point> points;

	for(int i=0;i<k;i++)
	{
		vector<double> point(dimensions,0);

		for(int j=0;j<dimensions;j++)
		{
			point[j]=rand()%MAX_VALUE;

		}
		clusters.push_back(Point(point));
	}

	x=(0.45/k)*MAX_VALUE;//min((0.001*sqrt(total_points)/(k*k)),

	int points_per_cluster=total_points/k;

	for(int i=0;i<k;i++)
	{
		vector<double> coords=clusters[i].get_coordinates();
		vector<double> temp(dimensions,0);
		for(int l =0; l<points_per_cluster; ++l)
		{
			for(int j=0;j<dimensions;j++)
			{
				double rand_01=((double)rand()/(double)RAND_MAX);


				int sign=rand()%2;
				if(sign==0)
				{
					temp[j] = rand_01*x+coords[j];
				}
				else
				{
					temp[j] = coords[j]-rand_01*x;
				}
			}
			points.push_back(Point(temp));
		}
	}
	return points;
}




/*
 * datatset.h
 *
 *  Created on: Nov 5, 2012
 *      Author: dishant
 */
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <math.h>

#include "kmeans.h"


#ifndef DATASET_H_
#define DATASET_H_

#define MIN_VALUE 0
#define MAX_VALUE 200000
/**
 * @arg1 dimensions of the input datatset
 * @arg2 clusters needed
 * @arg3 total number of points in the datatset > size of datatset
 */
vector<Point> get_datatset(int dimensions, int k, int total_points, vector<Cluster>& clusters, double& x)
{
	//vector<Cluster> clusters;
	vector<Point> points;

	for(int i=0;i<k;i++)
	{
		vector<double> point(dimensions,0);

		for(int j=0;j<dimensions;j++)
		{
			point[j]=rand()%MAX_VALUE;

		}
		clusters.push_back(Cluster(point));
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



#endif /* DATASET_H_ */


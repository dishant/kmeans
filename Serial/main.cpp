/*
 * main.cpp
 *
 *  Created on: Nov 5, 2012
 *      Author: dishant
 */


#include <iostream>
#include <vector>
#include "kmeans.h"
#include "dataset.h"


using namespace std;


int main()
{

	int dimensions = 3;
	int k = 100;
	int no_of_points = 200000;

	vector<Cluster> clusters;
	double x;

	vector<Cluster> c = kmeans(k,get_datatset(dimensions,k,no_of_points,clusters,x));
//	for(int i=0; i<c.size(); i++)
//	{
//		vector<double> d = c[i].get_coordinates();
//		cout<<"Output-Cluster"<<i<<" ("<<d[0]<<","<<d[1]<<')'<<endl;
//	}
//
//	for(int i=0; i<c.size(); i++)
//	{
//		vector<double> d = clusters[i].get_coordinates();
//		cout<<"Input-Cluster"<<i<<" ("<<d[0]<<","<<d[1]<<')'<<endl;
//	}

	cout<<endl<<x<<endl;

	cout<<check_accuracy(clusters,c,x)<<endl;

	return 0;
}




#include <iostream>
#include <vector>
#include "kmeans.h"
#include <mpi.h>

#define DIMENSION 3
#define NO_OF_POINTS 10000
#define K 20

using namespace std;

//end is actually end_index + 1.
double* serialize(vector<Point>& points, int start, int end, int& return_array_size)
{
		double num_points = end-start;
		return_array_size = 3+(num_points)*(DIMENSION+2);
		double* return_array = new double[return_array_size];
		int index = 0;


		return_array[index++] = 12345;
		return_array[index++] = (double) start;
		return_array[index++] = num_points;

		for(int i=0; i<num_points; ++i)
		{
			vector<double> coordinates = points[start+i].get_coordinates();
			for(int j=0; j<DIMENSION; ++j)
			{
				return_array[index++] = coordinates[j];
			}

			return_array[index++] = (double)points[start+i].get_cluster_id();
			return_array[index++] = points[start+i].get_cluster_distance();

		}

		return return_array;
}

vector<Point> deserialize(double* input, int& start_index)
{
		int index = 0;
		vector<Point> points;
		if(input[index++] != 12345) return points;
		start_index = input[index++];
		int num_points = input[index++];




		for(int i=0; i<num_points; ++i)
		{
				vector<double> coordinates;
				for(int j=0; j<DIMENSION; ++j)
				{
					coordinates.push_back(input[index++]);
				}

				points.push_back(Point(coordinates));
				points[i].set_cluster_id((int)input[index++]);
				points[i].set_cluster_distance(input[index++]);
		}

		return points;
}

int main(int argc, char **argv)
{

	MPI_Init (&argc, &argv);

	int mpi_id;//The id of this processor
	int mpi_size;// Total number of processors in MPI environment
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_id);

	vector<Point> points;//points contained by this processor are stored here
	vector<Point> clusters;//current cluster
	int start_index = 0; //indicates global index of points for this processor
	vector<Point> input_clusters;
	double x;

	if( mpi_id == 0 )
	{
		vector<Point> dataset;

		dataset = get_dataset(DIMENSION,K,NO_OF_POINTS,input_clusters,x);//NO_OF_POINTS may change after this step

		sort(dataset.begin(),dataset.end(),greater_than);
		int no_of_points = points.size();
		int points_per_cluster = no_of_points/K;


		//Create the clusters
		for(int j=0; j<points_per_cluster*K;j=j+points_per_cluster)
		{
			vector<double> sum_ppc(DIMENSION,0); //sum of points per cluster(ppc)
			for(int r=j; r<j+points_per_cluster;r++)
			{
				for(int z=0;z<DIMENSION;z++)
				{
					sum_ppc[z] += (dataset[r].get_coordinates())[z];
				}
			}

			for(int z=0;z<DIMENSION;z++)
			{
				sum_ppc[z] = sum_ppc[z]/points_per_cluster;
			}
			clusters.push_back(Point(sum_ppc));
		}

		//Communicate to others about the clusters
		double* send_array;
		int send_size;

		send_array = serialize(clusters,0,K,send_size);
		for(int i=1; i< mpi_size; ++i)
		{
			MPI_Send(send_array,send_size,MPI_DOUBLE,i,1,MPI_COMM_WORLD);
		}

		//Communicate to others about points
		int points_per_processor = no_of_points/mpi_size;

		for(int i=1; i<mpi_size; ++i)
		{

			if(i < (mpi_size -1))
				send_array = serialize(dataset,i*points_per_processor,(i+1)*points_per_processor,send_size);
			else
			{
				send_array = serialize(dataset,i*points_per_processor,no_of_points,send_size);
			}

			MPI_Send(send_array,send_size,MPI_DOUBLE,i,12345,MPI_COMM_WORLD);
		}

		/* Initialize processor-0 */
		for(int i=0; i<points_per_processor; ++i)
		{
			points.push_back(dataset[i]);
		}

	}

	else
	{
		int size_expected = ((NO_OF_POINTS/mpi_size)+mpi_size)*(DIMENSION+2)+3;
		double* receive_array = new double[size_expected];

		MPI_Status status;
		MPI_Recv(&receive_array,size_expected,MPI_DOUBLE,0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		clusters = deserialize(receive_array,start_index);//start-index is just 0

		MPI_Recv(&receive_array,size_expected,MPI_DOUBLE,0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		points = deserialize(receive_array,start_index);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	int t=0;

	while(t++<5)
	{
		//Associate the points to the clusters
		for(int j=0;j<points.size();j++)
		{
			  set_min_cluster(clusters,points[j]);
		}

		//Move the clusters according to the clusters and associated points
		vector<double> temp(DIMENSION,0);//just to create next variable

		/* Cluster-Id variable in point will be used to tell about no_of_points_in_cluster */
		vector<Point> sum_of_points_in_cluster(K,temp);


		for(int j=0;j<points.size();j++)
		{
			 (sum_of_points_in_cluster[points[j].get_cluster_id()]).add(points[j]);
			 int num_points_in_cluster = (sum_of_points_in_cluster[points[j].get_cluster_id()]).get_cluster_id();
			 (sum_of_points_in_cluster[points[j].get_cluster_id()]).set_cluster_id(num_points_in_cluster + 1);
		}

		double max_distance = 0;
		int farthest_point_index = -1;
		for(int j=0; j< points.size(); j++)
		{
			if(points[j].get_cluster_distance() > max_distance)
			{
				max_distance = points[j].get_cluster_distance();
				farthest_point_index = j;
			}
		}

		Point farthest_point(points[farthest_point_index]);
		farthest_point.set_cluster_distance(max_distance);

		/*
		 * Receive from other processors their sums and farthest points
		 * Send to all processors the new cluster points
		 */
		if(mpi_id == 0)
		{
			int size_expected = ((NO_OF_POINTS/mpi_size)+mpi_size)*(DIMENSION+2)+3;
			double* receive_array = new double[size_expected];

			MPI_Status status;

			for(int i=1; i<mpi_size; ++i)
			{
				MPI_Recv(&receive_array,size_expected,MPI_DOUBLE,i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				vector<Point> received_sum_of_points;
				received_sum_of_points = deserialize(receive_array,start_index);

				for(int j=0; j<K; ++j)
				{
					sum_of_points_in_cluster[j].add(received_sum_of_points[j]);
					int num_points_in_cluster = sum_of_points_in_cluster[j].get_cluster_id();
					sum_of_points_in_cluster[j].set_cluster_id(num_points_in_cluster + received_sum_of_points[j].get_cluster_id());
				}

				if(received_sum_of_points[K].get_cluster_distance()>farthest_point.get_cluster_distance())
				{
					farthest_point.set_coordinates(received_sum_of_points[K].get_coordinates());
					farthest_point.set_cluster_distance(received_sum_of_points[K].get_cluster_distance());
				}
			}

			bool farthest_point_assigned = false;
            for(int j=0;j<K;j++)
            {
            	int num_points_in_cluster = sum_of_points_in_cluster[j].get_cluster_id();
            	if(num_points_in_cluster > 0)
            	{
            		sum_of_points_in_cluster[j].divide_by(num_points_in_cluster);
					clusters[j].set_coordinates(sum_of_points_in_cluster[j].get_coordinates());
            	}
            	else if(!farthest_point_assigned)
            	{
            		clusters[j].set_coordinates(farthest_point.get_coordinates());
            		farthest_point_assigned = true;
            	}
            }

            //Communicate to others about the clusters
			double* send_array;
			int send_size;

			send_array = serialize(clusters,0,K,send_size);
			for(int i=1; i< mpi_size; ++i)
			{
				MPI_Send(send_array,send_size,MPI_DOUBLE,i,1,MPI_COMM_WORLD);
			}

		}

		/*
		 * Communicate to processor-0 about sums and farthest point
		 * Receive from processor-0 new cluster centers
		 */
		else
		{
			double* send_array;
			int send_size;
			//Piggy Back farthest point
			sum_of_points_in_cluster.push_back(farthest_point);
			send_array = serialize(sum_of_points_in_cluster,0,K+1,send_size);
			MPI_Send(send_array,send_size,MPI_DOUBLE,0,1,MPI_COMM_WORLD);

			int size_expected = ((NO_OF_POINTS/mpi_size)+mpi_size)*(DIMENSION+2)+3;
			double* receive_array = new double[size_expected];

			MPI_Status status;
			MPI_Recv(&receive_array,size_expected,MPI_DOUBLE,0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			clusters = deserialize(receive_array,start_index);//start-index is just 0
		}
		
	}

			if(mpi_id == 0)cout<<check_accuracy(clusters,input_clusters,x)<<endl;
	return 0;
}


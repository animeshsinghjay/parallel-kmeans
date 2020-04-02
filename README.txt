Directory structure:

	-openmp
		-make
		-./anim2
	-pthreads
		-make
		-./anim1
	-seq
		-make
		-./anim



Flags used during compilation:

-std=c++11
-fopenmp
-pthread




The main folder also has all the files.
To compile:

g++ main_omp.c lab1_io.c lab1_omp.cpp -std=c++11 -o anim2 -lpthread -fopenmp

To run:

./anim2 K T dataset_5000_4.txt outp_cluster2.txt outp_centroid2.txt



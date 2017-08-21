all: global_fit

global_fit: global_fit.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` -g -O3 -Wall -Wextra \
		global_fit.cc -o global_fit

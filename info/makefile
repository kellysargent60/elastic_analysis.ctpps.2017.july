all: print_info

print_info: print_info.cc
	g++ --std=c++11 `root-config --libs` `root-config --cflags` -g -O3 -Wall -Wextra \
		print_info.cc -o print_info

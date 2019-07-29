// Fingerprint.h
//
// Copyright (c) 2016
// Universitaet Duisburg-Essen
// Campus Duisburg
// Institut fuer Soziologie
// Prof. Dr. Rainer Schnell
// Lotharstr. 65
// 47057 Duisburg
//
// This file is part of the command line application "mbtSearch".
//
// "mbtSearch" is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// "mbtSearch" is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with "mbtSearch". If not, see <http://www.gnu.org/licenses/>.

#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>

// objects of class Timer provide simple functionalities
// to measure the elapsed time between two points in code

class Timer {
	private:
	
	struct timeval mStart;
	
	public:
	
	// start time
	inline void start() {
		gettimeofday(&mStart, NULL);
	}
	
	// return elpased time in ms since last start
	inline long stop() {
		struct timeval stop;
		
		gettimeofday(&stop, NULL);
		
		return (stop.tv_sec - mStart.tv_sec) * 1000 + (stop.tv_usec - mStart.tv_usec) / 1000;
	}
};
#endif

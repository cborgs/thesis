# Makefile
#
# Copyright (c) 2016
# Universitaet Duisburg-Essen
# Campus Duisburg
# Institut fuer Soziologie
# Prof. Dr. Rainer Schnell
# Lotharstr. 65
# 47057 Duisburg
#
# This file is part of the command line application "mbtSearch".
#
# "mbtSearch" is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# "mbtSearch" is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with "mbtSearch". If not, see <http://www.gnu.org/licenses/>.

CC=g++
CFLAGS=-c -O2 -Wall -pthread -Iinclude
LDFLAGS=-pthread
CPP_FILES = $(wildcard src/*.cpp)
OBJ_FILES = $(patsubst src/%.cpp,obj/%.o,$(CPP_FILES))
EXECUTABLE=mbtSearch

all: $(EXECUTABLE)

clean:
	rm -f $(EXECUTABLE) $(EXECUTABLE).exe $(OBJ_FILES)
	
$(EXECUTABLE): $(OBJ_FILES) 
	$(CC) $(LDFLAGS) $(OBJ_FILES) -o $@

obj/%.o: src/%.cpp include/*.h
	$(CC) $(CFLAGS) $< -o $@

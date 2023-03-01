
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wno-sign-compare -Ofast -march=native -fopenmp -g3
SRCS = main.cpp
OBJS = $(SRCS:.cpp=.o)
TARGET = microbude

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
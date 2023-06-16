#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <string>

template<typename units = std::chrono::nanoseconds>
struct timer {

	timer() {};

	void begin() {
		this->start_time = std::chrono::high_resolution_clock::now();
	}

	std::string display() {
		return std::to_string(std::chrono::duration_cast<units>(std::chrono::high_resolution_clock::now() - this->start_time).count());
	}

private:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;

};//timer

#endif// TIMER_HPP

/*
 * =============================================================================
 * stopwatch.h
 *
 * For timing information.
 * -----------------------------------------------------------------------------
 *
 * author : Marco Capettini <s252912@studenti.polito.it>
 *          Emanuele Parisi <emanuele.parisi@polito.it>
 *          Gianvito Urgese <gianvito.urgese@polito.it>
 * date   : December, 2019
 * =============================================================================
 */

#ifndef ISOMIR_SEA_STOPWATCH_H
#define ISOMIR_SEA_STOPWATCH_H

#include <chrono>

template<typename time_t = std::chrono::microseconds,
        typename clock_t=std::chrono::high_resolution_clock,
        typename duration_t=double>
class stopwatch
{
private:
    std::chrono::time_point<clock_t> _start, _end;
    std::string _label;
public:
    stopwatch(std::string str) {
        start();
        this->_label= str;
    }
    void start() {
        _start = _end = clock_t::now();
    }
    duration_t stop() {
        _end = clock_t::now(); return elapsed();
    }
    duration_t elapsed() {
        auto delta = std::chrono::duration_cast<time_t>(_end-_start);
        return delta.count();
    }
    void print(std::ofstream & out){
        this->stop();
        out << this->_label << " (Sec) = " << std::to_string(this->elapsed()/1000000) << std::endl;
    }
};

#endif //ISOMIR_SEA_STOPWATCH_H

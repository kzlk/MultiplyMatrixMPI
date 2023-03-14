#include <cstdlib>
#include <time.h>

typedef struct Timer {
    clock_t start;
    clock_t end;
} Timer;


Timer* createTimer() {
    Timer* timer = (Timer*)malloc(sizeof(Timer));
    timer->start = 0;
    timer->end = 0;
    return timer;
}

void resetTimer(Timer* timer) {
    timer->start = 0;
    timer->end = 0;
}

void startTimer(Timer* timer) {
    timer->start = clock();
}

void stopTimer(Timer* timer) {
    timer->end = clock();
}

double getElapsedSeconds(Timer* timer) {
    return ((double)(timer->end - timer->start)) / CLOCKS_PER_SEC;
}

void destroyTimer(Timer* timer) {
    free(timer);
}

#include "Timer.h"

timerify* createTimer() {
    timerify* timer = (timerify*)malloc(sizeof(timer));
    timer->start = 0;
    timer->end = 0;
    return timer;
}

void resetTimer(timerify* timer) {
    timer->start = 0;
    timer->end = 0;
}

void startTimer(timerify* timer) {
    timer->start = clock();
}

void stopTimer(timerify* timer) {
    timer->end = clock();
}

double getElapsedSeconds(timerify* timer) {
    return ((double)(timer->end - timer->start)) / CLOCKS_PER_SEC;
}

void destroyTimer(timerify* timer) {
    free(timer);
}

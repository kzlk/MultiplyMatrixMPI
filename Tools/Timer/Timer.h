#pragma once

#include <cstdlib>
#include <time.h>

typedef struct Timer {
    clock_t start;
    clock_t end;
} timerify;

timerify* createTimer();
void resetTimer(timerify* timer);
void startTimer(timerify* timer);
void stopTimer(timerify* timer);
double getElapsedSeconds(timerify* timer);
void destroyTimer(timerify* timer);

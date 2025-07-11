#ifndef _BREATHPARTICLES_H_
#define _BREATHPARTICLES_H_

#include <vector>

#include "variables.h"

void breathAddParticles(double time, std::shared_ptr< UserCtx > & user,
                        std::vector<Particle> & Ptc);

void breathCalcPartDia(double time, Particle & P);

void breathCalcPartDia(double time, std::vector<Particle> & Ptc);

#endif

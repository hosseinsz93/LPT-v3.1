#ifndef _REYTOCD_H
#define _REYTOCD_H

#include <array>
#include <cmath>


template <typename T>
T interpolate(T x0, T x, T x1, T y0, T y1)
{
  T ratio = (y1 - y0) / (x1 - x0);
  return y0 + ratio * (x - x0);
}

class ReysToCd
{
  private:
    static constexpr std::array<double, 44> reys {
                        { 0.111033632, 0.311534748, 1.251418576, 3.258269574,
                          10.15064189, 23.44929931, 71.96856730, 224.2070945,
                          510.2608026, 861.1224963, 1431.674059, 2111.909602,
                          2890.939220, 8233.498592, 17130.34990, 27641.29854,
                          49522.80037, 87409.46087, 117876.8635, 173884.0473,
                          217601.7269, 248945.2798, 264288.2907, 280576.9229,
                          306911.3765, 330735.2764, 345908.3446, 395733.3007,
                          466478.0026, 533669.9231, 629073.3667, 730527.1543,
                          848342.8982, 1061632.062, 1348559.553, 1738840.473,
                          2208797.205, 2848035.868, 3357175.352, 4393970.561,
                          5581529.787, 6881170.798, 7990931.404, 9561355.820 }};

    static constexpr std::array<double, 44> cd {
                        { 259.8259595, 86.79447159, 24.15107730, 10.13800220,
                          4.454608319, 2.459648318, 1.268161613, 0.716405201,
                          0.557226480, 0.485856410, 0.433415823, 0.404708995,
                          0.369368893, 0.361027000, 0.377903533, 0.395568975,
                          0.414060204, 0.453676237, 0.443430332, 0.423627483,
                          0.369368893, 0.329501308, 0.256289113, 0.208662518,
                          0.155051578, 0.120600527, 0.098189148, 0.089615050,
                          0.087591168, 0.089615050, 0.110069417, 0.115214711,
                          0.132139222, 0.155051578, 0.173811843, 0.199344003,
                          0.228626719, 0.239314082, 0.250501035, 0.274468219,
                          0.287298483, 0.287298483, 0.293936804, 0.300728510 }};

  public:
    static double conv(double reysArg)
    {
      if (reysArg <= 1.0e-15)
        return 0.0;

      if (reysArg <= 0.1)
        return 24.0 / reysArg;

      if (reysArg <= reys[1])
      {
        auto r = interpolate(log(reys[0]), log(reysArg), log(reys[1]),
                             log(cd[0]),                 log(cd[1]));
        return exp(r);
      }

      unsigned ih = std::lower_bound(reys.begin(), reys.end(), reysArg)
                                                                - reys.begin();
      if (ih >= reys.size())
        ih = reys.size() - 1;
      int il = ih - 1;

      auto r = interpolate(log(reys[il]), log(reysArg), log(reys[ih]),
                           log(cd[il]),                 log(cd[ih]));
      return exp(r);
    }
};

#endif

#ifndef FIRE_H
#define FIRE_H


class Fire : public Minimiser {
  public:
    Fire(State &state);
    Fire(State &state, AdjustFunc adjustModel);
    ~Fire() {};

    Fire& setDtMax(double dt_max);

    void iteration();

    override bool checkConvergence();

  private:
    int _n_min = 5;
    double _f_inc = 1.1;
    double _f_dec = 0.5;
    double _f_a = 0.99;
    double _a_start = 0.1;
    double _a;
    double _dt;
    double _dt_max;
    double _gnorm;
    std::vector<double> _g;
    std::vector<double> _v;

    void getDirection();
};

#endif

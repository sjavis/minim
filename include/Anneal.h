#ifndef ANNEAL_H
#define ANNEAL_H


class Anneal : public Minimiser {
  public:
    Anneal(State &state, double temp_init, double displacement, AdjustFunc adjustModel=NULL);
    ~Anneal() {};

    Anneal& setMaxIter(int maxIter);
    Anneal& setTempInit(double temp_init);
    Anneal& setDisplacement(double displacement);

    void iteration();
    bool checkConvergence() override;

  private:
    int _since_accepted;
    int _max_rejections;
    double _temp;
    double _temp_init;
    double _displacement;
    double _current_e;
    std::vector<double> _current_state;

    bool acceptMetropolis(double energy);
};

#endif

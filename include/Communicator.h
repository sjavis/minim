#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>

class Args;
class Priv;


class Communicator {
  public:
    Communicator(int ndof, Args &args);
    ~Communicator();

    std::vector<double> assignBlock(const std::vector<double> &in);

    void communicate(std::vector<double> &vector);
    double get(const std::vector<double> &vector, int i);
    double dotProduct(const std::vector<double> &a, const std::vector<double> &b);

    std::vector<double> gather(const std::vector<double> &block, int root=-1);
    std::vector<double> scatter(const std::vector<double> &data, int root=-1);
    void bcast(double &value, int root=0);
    void bcast(std::vector<double> &value, int root=0);

  private:
    Priv *priv;
};

#endif

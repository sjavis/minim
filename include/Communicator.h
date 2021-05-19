#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

#include <vector>

class Args;
class Priv;


class Communicator {
  public:
    Communicator(int ndof, Args &args);
    ~Communicator();

    std::vector<double> assignBlock(std::vector<double> in);

    void communicate(std::vector<double> &vector);
    double get(std::vector<double> vector, int i);
    double dotProduct(std::vector<double> a, std::vector<double> b);
    std::vector<double> gather(std::vector<double> block, int ndof);

  private:
    Priv *priv;
};

#endif

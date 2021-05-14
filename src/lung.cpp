#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <memory>
#include <fcntl.h>
#include <inttypes.h>

class Random {
 private:
  std::mt19937_64 generator;

 public:
  Random(unsigned seed)
      : generator(seed) {}

  int get(int64_t begin, int64_t end) {
    return std::uniform_int_distribution<int64_t>(begin, end - 1)(generator);
  }

};

struct Level {
  int L = 0, r = 0, count = 0;
  double bAngle = 0.0, gAngle = 0.0;
  int centroid[3] = {0, 0, 0};
  double direction[3] = {0.0, 0.0, 0.0};
};

int64_t numAirways = 0, numAirwayCells = 0, numOutOfBoundsCells = 0;
std::set<int64_t> airwayEpiCellPositions1D;
std::vector<Level> levels;
std::shared_ptr<Random> rnd_gen = std::make_shared<Random>(753695190);
#define DIMENSION (300)
int gridSize[3] = {DIMENSION, DIMENSION, DIMENSION};

void rotate(int (&vec)[3], const double (&axis)[3], double angle) {
  /**
   * http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html
   */
  double cosAngle = cos(angle);
  double sinAngle = sin(angle);
  double x = vec[0], y = vec[1], z = vec[2];
  // double u = 0, v = 1, w = 0;
  double u = axis[0], v = axis[1], w = axis[2];
  double dotProduct = (x * u) + (y * v) + (z * w);
  double newX = u * dotProduct * (1 - cosAngle) + x * cosAngle + (-w * y + v * z) * sinAngle;
  double newY = v * dotProduct * (1 - cosAngle) + y * cosAngle + (w * x - u * z) * sinAngle;
  double newZ = w * dotProduct * (1 - cosAngle) + z * cosAngle + (-v * x + u * y) * sinAngle;
  vec[0] = (int)round(newX);
  vec[1] = (int)round(newY);
  vec[2] = (int)round(newZ);
}

void constructSegment(const int (&root)[3],
  const Level &level,
  double rotateZ,
  bool isTerminal,
  int (&newRoot)[3]) {
  // Build cylinder at origin along y-axis
  double cylinderRadialIncrement = M_PI / 2;
  for (int z = 0; z <= level.L; z++) {
    for (double az = 0; az < 2 * M_PI; az += M_PI / 180) {
#ifndef COMPUTE_ONLY
      int x = (int)round(level.r * sin(cylinderRadialIncrement) * cos(az));
      int y = (int)round(level.r * sin(cylinderRadialIncrement) * sin(az));
      int position[3] = {x, y, z};
      // Treat as positional vectors and apply rotations and translate y
      rotate(position, {0.0, 1.0, 0.0}, level.bAngle);
      // Rotate z
      rotate(position, {0.0, 0.0, 1.0}, rotateZ);
      position[0] += root[0];
      position[1] += root[1];
      position[2] += root[2];
      // Set and verify new location is within grid
      if ((0 <= position[0] && position[0] < gridSize[0]) &&
          (0 <= position[1] && position[1] < gridSize[1]) &&
          (0 <= position[2] && position[2] < gridSize[2])) {
            airwayEpiCellPositions1D.insert(position[0] + position[1]
              * gridSize[0] + position[2] * gridSize[0] * gridSize[1]);
            numAirwayCells++;
      } else {
        numOutOfBoundsCells++;
      }
#else
      numAirwayCells++;
#endif
    }
  }
  // Create root for next generation
  int base[3] = {0, 0, level.L};
  rotate(base, {0.0, 1.0, 0.0}, level.bAngle);
  rotate(base, {0.0, 0.0, 1.0}, rotateZ);
  newRoot[0] = root[0] + base[0];
  newRoot[1] = root[1] + base[1];
  newRoot[2] = root[2] + base[2];
  numAirways++;
}

void construct(int (&root)[3],
  int iteration,
  int index,
  int end,
  double previousBranchAngle,
  double previousRotAngle) {
  if (iteration > end) {
    return;
  }
  bool isTerminal = (iteration == end) ? true : false;
  // Uniform randomly rotate each branch
  double rotateZ = (iteration >= 2) ? rnd_gen->get(0, 180) * M_PI / 180 : 0.0;
  rotateZ = previousRotAngle - rotateZ;
  // Draw left child bronchiole
  Level lvl = levels.at(index);
  lvl.bAngle = previousBranchAngle - lvl.bAngle;
  int child[3] = {0, 0, 0};
  constructSegment(root, lvl, rotateZ, isTerminal, child);
  construct(child, iteration + 1, index + 1, end, lvl.bAngle, rotateZ);
  // Uniform randomly rotate each branch
  rotateZ = (iteration >= 2) ? rnd_gen->get(0, 180) * M_PI / 180 : 0.0;
  rotateZ = previousRotAngle + rotateZ;
  // Draw right child bronchiole
  lvl = levels.at(index);
  lvl.bAngle = previousBranchAngle + lvl.bAngle;
  lvl.gAngle = -lvl.gAngle;
  constructSegment(root, lvl, rotateZ, isTerminal, child);
  construct(child, iteration + 1, index + 1, end, lvl.bAngle, rotateZ);
}

void loadEstimatedParameters(std::vector<Level>& levels) {
  /* Yeh et al 1980
    scale = 2000  => 1 unit = 5um, 1cm = 10^-2m = 10^4 um, 10^4/5 = 2000 units
  */
  double scale = 2000;
  std::ifstream file;
  file.open("table.txt");
  std::string line;
  if (file.is_open()) {
    while (std::getline(file, line)) {
      std::stringstream lstream(line);
      Level e;
      int tempI;
      double tempD;
      lstream >> tempI;
      lstream >> e.count;
      lstream >> tempD;
      e.L = (int)round(scale * tempD);
      lstream >> tempD;
      e.r = (int)(round(scale * tempD) / 2);
      lstream >> tempI;
      e.bAngle = tempI * M_PI / 180;
      lstream >> tempI;
      e.gAngle = tempI * M_PI / 180;
      if (levels.size() > 73) {  // Negate for left lobe data
        e.bAngle *= -1.0;
      }
      levels.push_back(e);
    }
    std::printf("Loaded estimated levels %ld\n", levels.size());
    file.close();
  } else {
    std::printf("Failed to open table.txt\n");
  }
}

int main(int argc, char *argv[]) {
  loadEstimatedParameters(levels);
  // Generate partial upper right lobe starting at root and recursively build tree
  int generations = 3;
  int startIndex = 24 - generations;
  Level lvl = levels.at(startIndex);
  int root[3] = {gridSize[0]/6, gridSize[1]/2, 0};
  int child[3] = {0, 0, 0};
  constructSegment(root, lvl, 0.0, false, child);
  construct(child, 1, startIndex + 1, generations - 1, lvl.bAngle, 0.0);
  // Print stats
  std::printf("Airways " "%" PRId64 "\n", numAirways);
  std::printf("Epithelial cells generated " "%" PRId64, numAirwayCells);
  std::printf(" Out of bounds " "%" PRId64 "\n", numOutOfBoundsCells);
  // Write airway epithileal cells
  std::ofstream ofs;
  ofs.open("airway.csv", std::ofstream::out | std::ofstream::app);
  if (!ofs.is_open()) {
    std::printf("Could not create file");
    exit(1);
  }
  for (const int64_t& a : airwayEpiCellPositions1D) {
    ofs << a << std::endl;
  }
  std::printf("Written %ld", (long)ofs.tellp());
  std::printf(" Stored " "%" PRId64 "\n", airwayEpiCellPositions1D.size());
  ofs.close();
  return 0;
}

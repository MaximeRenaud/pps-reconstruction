#ifndef GastofCommons_h
#define GastofCommons_h

#include "TString.h"
#include <iostream>
#include <map>

class Coord {
  public:
    inline Coord() : x(0), y(0) {;}
    inline Coord(unsigned int x_, unsigned int y_) : x(x_), y(y_) {;}
    inline ~Coord() {;}

    bool operator==(const Coord& c) const { return (x==c.x && y==c.y); }
    //friend std::ostream& operator<<(std::ostream& os, const Coord& c);
    inline void Print() { std::cout << Form("(%2d, %2d)", x, y) << std::endl; }
    
    unsigned int x;
    unsigned int y;
};

/*std::ostream&
operator<<(std::ostream& os, const Coord& c)
{
  os << "(" << c.x << ", " << c.y << ")";
  return os;
  }*/

class GastofCoordinatesMap {
  public:
    inline ~GastofCoordinatesMap() { fBuilt = false; fMap.clear(); }
    inline static Coord GetCoordinates(unsigned int i) {
      if (i>63) return Coord();
      return GetInstance()->fMap[i];
    }
    inline static bool IsNeighbour(const Coord& c1, const Coord& c2) {
      if (c1==c2) return false;
      //if (abs(c1.x-c2.x)<2 && abs(c1.x-c2.x)<2) return true;
      if (c1.x<8 && c2.x==c1.x+1 && c2.y==c1.y) return true;
      if (c1.x>1 && c2.x==c1.x-1 && c2.y==c1.y) return true;
      if (c1.y<8 && c2.x==c1.x && c2.y==c1.y+1) return true;
      if (c1.y>1 && c2.x==c1.x && c2.y==c1.y-1) return true;
      return false;
    }
    inline static bool IsNeighbour(unsigned int id1, unsigned id2) {
      return IsNeighbour(GetCoordinates(id1), GetCoordinates(id2));
    }
    inline static int GetChannelId(const Coord& c) {
      map<unsigned int, Coord> m = GetInstance()->fMap;
      for (map<unsigned int, Coord>::const_iterator it=m.begin(); it!=m.end(); it++) {
        if (it->second==c) return it->first;
      }
      return -1;
    }
    typedef map<unsigned int, Coord> Neighbours;
    inline static Neighbours GetNeighbours(const Coord& c1) {
      Neighbours out;
      if (c1.x<8) { Coord c(c1.x+1, c1.y); out.insert(pair<unsigned int, Coord>(GetChannelId(c), c)); }
      if (c1.x>1) { Coord c(c1.x-1, c1.y); out.insert(pair<unsigned int, Coord>(GetChannelId(c), c)); }
      if (c1.y<8) { Coord c(c1.x, c1.y+1); out.insert(pair<unsigned int, Coord>(GetChannelId(c), c)); }
      if (c1.y>1) { Coord c(c1.x, c1.y-1); out.insert(pair<unsigned int, Coord>(GetChannelId(c), c)); }
      return out;
    }
    inline static GastofCoordinatesMap* GetInstance() {
      if (fBuilt) return fCM;
      fCM = new GastofCoordinatesMap;
      return fCM;
    }
  private:
    static GastofCoordinatesMap* fCM;
    inline GastofCoordinatesMap() {
      fMap[24] = Coord(1, 8); fMap[25] = Coord(2, 8); fMap[26] = Coord(3, 8); fMap[27] = Coord(4, 8); fMap[28] = Coord(5, 8); fMap[29] = Coord(6, 8); fMap[30] = Coord(7, 8); fMap[31] = Coord(8, 8);
      fMap[20] = Coord(1, 7); fMap[21] = Coord(2, 7); fMap[60] = Coord(3, 7); fMap[61] = Coord(4, 7); fMap[62] = Coord(5, 7); fMap[63] = Coord(6, 7); fMap[22] = Coord(7, 7); fMap[23] = Coord(8, 7);
      fMap[18] = Coord(1, 6); fMap[54] = Coord(2, 6); fMap[55] = Coord(3, 6); fMap[56] = Coord(4, 6); fMap[57] = Coord(5, 6); fMap[58] = Coord(6, 6); fMap[59] = Coord(7, 6); fMap[19] = Coord(8, 6);
      fMap[16] = Coord(1, 5); fMap[48] = Coord(2, 5); fMap[49] = Coord(3, 5); fMap[50] = Coord(4, 5); fMap[51] = Coord(5, 5); fMap[52] = Coord(6, 5); fMap[53] = Coord(7, 5); fMap[17] = Coord(8, 5);
      fMap[14] = Coord(1, 4); fMap[42] = Coord(2, 4); fMap[43] = Coord(3, 4); fMap[44] = Coord(4, 4); fMap[45] = Coord(5, 4); fMap[46] = Coord(6, 4); fMap[47] = Coord(7, 4); fMap[15] = Coord(8, 4);
      fMap[12] = Coord(1, 3); fMap[36] = Coord(2, 3); fMap[37] = Coord(3, 3); fMap[38] = Coord(4, 3); fMap[39] = Coord(5, 3); fMap[40] = Coord(6, 3); fMap[41] = Coord(7, 3); fMap[13] = Coord(8, 3);
      fMap[ 8] = Coord(1, 2); fMap[ 9] = Coord(2, 2); fMap[32] = Coord(3, 2); fMap[33] = Coord(4, 2); fMap[34] = Coord(5, 2); fMap[35] = Coord(6, 2); fMap[10] = Coord(7, 2); fMap[11] = Coord(8, 2);
      fMap[ 0] = Coord(1, 1); fMap[ 1] = Coord(2, 1); fMap[ 2] = Coord(3, 1); fMap[ 3] = Coord(4, 1); fMap[ 4] = Coord(5, 1); fMap[ 5] = Coord(6, 1); fMap[ 6] = Coord(7, 1); fMap[ 7] = Coord(8, 1);
      fBuilt = true;
    }
    static bool fBuilt;
    map<unsigned int, Coord> fMap;
};

bool GastofCoordinatesMap::fBuilt = false;
GastofCoordinatesMap* GastofCoordinatesMap::fCM = 0;

inline Coord GetCoordinates(unsigned short nino_id, unsigned short channel_id) {
  return GastofCoordinatesMap::GetCoordinates(nino_id*32+channel_id);
}

#endif

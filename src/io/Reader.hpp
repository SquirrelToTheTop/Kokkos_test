#include <string>
#include <cassert>

class Reader;

class Reader {

  public:
    Reader(){};
    ~Reader(){};

    // initialize Reader
    virtual bool initializeReader() = 0;

    // open object
    virtual bool open( const std::string &dir ) = 0;

  private:
    
};
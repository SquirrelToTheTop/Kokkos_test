#include "Reader.hpp"
#include "../datamodel/lightAMR.hpp"

// Hercule
#include "HIc.h"

HIC_USE;

class Reader_Hercule;

class Reader_Hercule : public Reader{

  public:
    Reader_Hercule(){
      _apiInitialized = false;
    };

    ~Reader_Hercule() = default;

    /**
     * Initialize API
    */
    bool initializeReader();

    /**
     * Open files or database (data container) 
     */
    bool open( const std::string &dir );

    LightAMR getAMRData( int tid, int did, const std::string &objName );

  private:

    /* Private variables for this class ---------------------------------------------- */

    // Api et base
    HIc_Api _hic_api;
    HIc_Base _hic_base;

    // ctx qui changera au fur et a mesure des appels pour traiter les domaines
    HIc_Ctx _hic_ctx;

    std::map<std::string, HIc_Obj> _loaded_object;

    // liste des indices des temps (en temps de la simulation en floats)
    std::vector<double> _base_times_indexes;

    bool _apiInitialized;

    /* Private functions for this class ---------------------------------------------- */

    /**
     * Permet d'acceder a une variable d'un objet XML de Hercule
    */
    template<typename var_t>
    inline
    var_t _hic_getVar( HIc_Obj &hicObj, const std::string &varName ) const {
      var_t tmp = static_cast<var_t>( 0 );
      try{
          hicObj.getAttrVal( varName, tmp );
      } catch ( std::exception &_exp ){
          std::cerr << "\t> [accessVar] Failed for '" << varName << "'" << std::endl;
          std::cerr << "\t\t> Error: " << _exp.what() << std::endl;
      }

      return tmp;
    }

};
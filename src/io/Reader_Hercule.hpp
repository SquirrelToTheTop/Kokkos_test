#include "Reader.hpp"
#include "../datamodel/lightAMR.hpp"

// Hercule
#include "HIc.h"

// for now
#define DEBUG_READER

HIC_USE;

class Reader_Hercule;

class Reader_Hercule : public Reader{

  public:

    Reader_Hercule(){
      _apiInitialized = false;
      _base_ndomains  = 0;
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

    /**
     * Retourne le nombre de sous-domaine contenu dans la base
    */
    inline uint32_t GetNumberOfDomains() const{ return _base_ndomains; }

    /**
     * tid (in): time-index of time to load
     * did (in): domain-index to load within the time
    */
    LightAMR GetAMRData( int tid, int did, const std::string &objName, int dim );

  private:

    /* Private variables for this class ---------------------------------------------- */

    // Api et base
    HIc_Api _hic_api;
    HIc_Base _hic_base;

    // ctx qui changera au fur et a mesure des appels pour traiter les domaines
    HIc_Ctx _hic_ctx;

    std::map<std::string, HIc_Obj> _loaded_object;

    // liste des indices des temps (en temps de la simulation en floats)
    uint32_t _base_ndomains;
    std::vector<double> _base_times_indexes;

    bool _apiInitialized;

    /* Private functions for this class ---------------------------------------------- */

    /**
     * Permet d'acceder a un attribut d'un objet XML de Hercule.
     * 
     * @param hicObj   (in): objet hercule portant l'attribut
     * @param attrName (in): identifiant (unique) de l'attribut
    */
    template<typename var_t>
    inline
    var_t _hic_getVar( HIc_Obj &hicObj, const std::string &attrName ) const {

      var_t tmp = static_cast<var_t>( 0 );
      
      try{
          hicObj.getAttrVal( attrName, tmp );
      } catch ( std::exception &_exp ){
          std::cerr << "\t> [Reader_Hercule::_hic_getVar] Failed for '" << attrName << "'" << std::endl;
          std::cerr << "\t\t> Error: " << _exp.what() << std::endl;
      }

      return tmp;
    }

    /**
     * Permet d'acceder a un attribut de type tableau d'un objet XML de Hercule en remplissant un tableau
     * pré-alloué.
     * 
     * @param hicObj   (in): objet hercule portant l'attribut
     * @param attrName (in): identifiant (unique) de l'attribut
     * @param nelem    (in): nombre d'élements de l'attribut
    */
    template<typename var_t>
    inline
    void _hic_getVarTab1d( HIc_Obj &hicObj, const std::string &attrName, DataArray1D_host<var_t> &arr ) const {
      assert( arr.extent(0) > 0 );

      try{
          hicObj.getAttrVal( attrName, arr.data(), arr.size() ); // nElemRead );
      } catch ( std::exception &_exp ){
          std::cerr << "\t> [accessVar] Failed for '" << attrName << "'" << std::endl;
          std::cerr << "\t\t> Error: " << _exp.what() << std::endl;  
      }
    }



};
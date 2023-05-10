#include "Reader_Hercule.hpp"

bool Reader_Hercule::initializeReader(){
  bool success = true;

  // check directory exists

  _hic_api = HIc_Api("API", "stdCommon");
  HIc_Init_Services( _hic_api );                   // usefull ? 
  HIc_Init_Standard_Services( _hic_api );          // usefull ? 
  HIc_Init_Standard_Site( _hic_api, "stdCommon" ); // usefull ? 

  if( _hic_api.isNull() ) success = false;

  if( success ) _apiInitialized = true;

  return success;
}

/**
 * 
 * @param bdir (in): chemin vers le dossier content la base HDep 
*/
bool Reader_Hercule::open( const std::string &bdir ) {

  bool success = true;

  if( ! _apiInitialized ){
    std::cerr << "\t> [Reader_Hercule]::open {FAILURE} API not initialized " << std::endl;
    return false;
  }

  _hic_base = HIc_Base( _hic_api, "depouillement", "parallele" );

  // configuration de la base Hercule
  _hic_base.setItemConf( "mode", "read" );
  _hic_base.setItemConf( "read_dir", bdir );
  _hic_base.setItemConf( "bd_name", "HDep" );
  _hic_base.setItemConf( "access_mode", "sequential_access" );
  _hic_base.open();

  // c'est necessaire avec Hercule pour les acces apres
  _hic_base.getTimeList( _base_times_indexes ); 

  return _hic_base.isOpen();
}

LightAMR Reader_Hercule::getAMRData( int tit, int did, const std::string &objName ){

  LightAMR lamr;

  if( _hic_ctx.isNull() ){
    // assert( tit < static_cast<int>( _base_times_indexes.size() ) );
    _hic_ctx = _hic_base.getCtxPar( _base_times_indexes[ tit ], did );
    _hic_ctx.open();
  }

  HIc_Obj ob_root = _hic_ctx.getRoot();
  HIc_Obj ob = ob_root.searchUniq( objName );



  return std::move( lamr );
}
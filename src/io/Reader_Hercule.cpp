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

  // On garde une trace du nombre de pas de temps + du nombre de domaines
  _hic_base.getTimeList( _base_times_indexes ); 
  _base_ndomains = _hic_base.getNbDomains();

#ifdef DEBUG_READER
  std::cout << "\t> [INFO][Reader_Hercule] Base contains " << _base_times_indexes.size() << " contextes \n";
  std::cout << "\t> [INFO][Reader_Hercule] Base contains " << _base_ndomains << " domains \n";
#endif

  return _hic_base.isOpen();
}

LightAMR Reader_Hercule::GetAMRData( int tit, int did, const std::string &objName ) {

  LightAMR lamr;

  if( _hic_ctx.isNull() ){
    // assert( tit < static_cast<int>( _base_times_indexes.size() ) );
    _hic_ctx = _hic_base.getCtxPar( _base_times_indexes[ tit ], did );
    _hic_ctx.open();
  }

  // Récupère l'objet d'information sur la simu
  HIc_Obj ramsesInfo = _hic_ctx.searchUniq( "RamsesInfo" );
  bool compressed = static_cast<bool>( _hic_getVar<int>( ramsesInfo, "isCompressed" ) );
  lamr.setAMRCompressionState( compressed );

  HIc_Obj ob_root = _hic_ctx.getRoot();
  HIc_Obj ob_amr = ob_root.searchUniq( objName );

  // nombre de cellules (modèle lightAMR)
  int64_t nbCells = _hic_getVar<int64_t>( ob_amr, "nbElements" );
  lamr.setNumberOfCells( nbCells );

  // tableau de raffinement & d'ownership
  auto & refArr = lamr.getRefinementArrayPtr();
  auto & ownArr = lamr.getOwnershipArrayPtr();

  // en CPS52 la taille des tableaux est != nbElements et 
  // n'est pas forcément la même pour les deux
  if( compressed ) {

    int64_t nbChars = _hic_getVar<int64_t>( ob_amr, "nbParentChars" );
    Kokkos::resize( refArr, nbChars );
    _hic_getVarTab1d<uint8_t>( ob_amr, "isParentChars", refArr );

    nbChars = _hic_getVar<int64_t>( ob_amr, "nbMaskChars" );
    Kokkos::resize( ownArr, nbChars );
    _hic_getVarTab1d<uint8_t>( ob_amr, "isMaskChars", ownArr );

  }else{

    Kokkos::resize( refArr, nbCells );
    _hic_getVarTab1d<uint8_t>( ob_amr, "isParentInt", refArr );

    Kokkos::resize( ownArr, nbCells );
    _hic_getVarTab1d<uint8_t>( ob_amr, "isMaskInt", ownArr );

  }

  uint64_t nlevels = _hic_getVar<uint64_t>( ob_amr, "nbLevels" );

  auto & ncplArr = lamr.getNcplArrayPtr();
  Kokkos::resize( ncplArr, nlevels );
  _hic_getVarTab1d<uint64_t>( ob_amr, "nbElementsPerLevel", ncplArr );

#ifdef DEBUG_READER
  std::cout << "\t> [INFO][Reader_Hercule]: load 'nbElements' {" << ob_amr.getAttrTypeName( "nbElements" ) 
            << "} to {" << getCustomTypeName( nbCells ) << "}, value : " << nbCells << std::endl;
  
  if( compressed ){
    std::cout << "\t> [INFO][Reader_Hercule]: load 'isParentChars' {" << ob_amr.getAttrTypeName( "isParentChars" ) 
              << "} to {" << getCustomTypeName( refArr ) << "}" << std::endl;
    std::cout << "\t> [INFO][Reader_Hercule]: load 'isMaskChars' {" << ob_amr.getAttrTypeName( "isMaskChars" ) 
              << "} to {" << getCustomTypeName( ownArr ) << "}" << std::endl;
  }else{
    std::cout << "\t> [INFO][Reader_Hercule]: load 'isParentInt' {" << ob_amr.getAttrTypeName( "isParentInt" ) 
              << "} to {" << getCustomTypeName( refArr ) << "}" << std::endl;
    std::cout << "\t> [INFO][Reader_Hercule]: load 'isMaskInt' {" << ob_amr.getAttrTypeName( "isMaskInt" ) 
              << "} to {" << getCustomTypeName( ownArr ) << "}" << std::endl;
  }

  std::cout << "\t> [INFO][Reader_Hercule]: load 'nbLevels' {" << ob_amr.getAttrTypeName( "nbLevels" ) 
            << "} to {" << getCustomTypeName( nlevels ) << "}, value : " << nlevels << std::endl;
    std::cout << "\t> [INFO][Reader_Hercule]: load 'nbElementsPerLevel' {" << ob_amr.getAttrTypeName( "nbElementsPerLevel" ) 
            << "} to {" << getCustomTypeName( ncplArr ) << "}" << std::endl;
#endif

  return lamr;
}
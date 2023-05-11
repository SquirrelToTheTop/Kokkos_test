#include "../utils/hzlnt.hpp"
#include "Kokkos_UnorderedMap.hpp"

#include "lightAMR.hpp"

#include <unordered_map>

class Collection;

using Collection_container = Kokkos::UnorderedMap<uint32_t, LightAMR>;
// using Collection_container = std::unordered_map<uint32_t, LightAMR>;

class Collection{
  public:

    Collection(){};
    Collection( const Collection & ) = delete;

    ~Collection(){};

    /**
     * @param key (in): identifiant unique pour l'élement LightAMR dans la collection (ID du domaine)
     * @param rv  (in): objet de type LightAMR
    */
    inline void addItem( uint32_t key, LightAMR && rv ) noexcept {
      if( _collections.exists( key ) ){
        std::cerr << "\t>[Collection]::addItem key ('" << key << "') already exists !" << std::endl;
      }else{
        _collections.insert( key, rv );
      }
    }

    inline uint32_t GetNumberOfItems() const { return static_cast<uint32_t>( _collections.size() ); }
  private:
    Collection_container _collections;
  
};
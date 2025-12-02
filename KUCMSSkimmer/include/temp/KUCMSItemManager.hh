// -*- C++ -*-
//
//  ItemManager Template Implementation
//  Jack W King III
//  Created:  Wed, 27 Jan 2021
//

#ifndef KUCMSItemManagerHeader_hh
#define KUCMSItemManagerHeader_hh

// basic C++ types
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <iostream>
#include <limits>

// Root includes
#include "TTree.h"

// KUCMS includes
#include "KUCMSItem.hh"

//......................................................................................
//   ItemManager Template Class
//......................................................................................

template < class T , template< class > class C = Item >
class ItemManager {

public:

    void set( std::string key, std::string name, std::string doc );
    void set( std::string name );
    void set( std::string name, T val );

    void fill( std::string key, T val );
    T get( std::string key );

    void clear();
    void clear( std::string key );
    void reset();

    T operator()( std::string key );

private:

    std::map< std::string, C<T> > items;
    bool valid( std::string key );
};

//////////////////////////////////////////////////////////////////////////////////////////
// Implementation must be in header (template)
//////////////////////////////////////////////////////////////////////////////////////////

template < class T , template< class > class C >
void ItemManager<T,C>::set( std::string key, std::string name, std::string doc )
{
    C<T> newitem( name, doc );
    items[key] = newitem;
}

template < class T , template< class > class C >
void ItemManager<T,C>::set( std::string name )
{
    set( name, name, "" );
}

template < class T , template< class > class C >
void ItemManager<T,C>::set( std::string name, T val )
{
    C<T> newitem( name, val, "" );
    items[name] = newitem;
}

template < class T , template< class > class C >
void ItemManager<T,C>::fill( std::string key, T val )
{
    if( valid( key ) ) items[key].fill( val );
}

template < class T , template< class > class C >
T ItemManager<T,C>::get( std::string key )
{
    return valid( key ) ? items[key].getvalue() : T();
}

template < class T , template< class > class C >
void ItemManager<T,C>::clear()
{
    for( auto & item : items ) item.second.clear();
}

template < class T , template< class > class C >
void ItemManager<T,C>::reset()
{
    for( auto & item : items ) item.second.erase();
}

template < class T , template< class > class C >
void ItemManager<T,C>::clear( std::string key )
{
    if( valid( key ) ) items[key].clear();
}

template < class T , template< class > class C >
T ItemManager<T,C>::operator()( std::string key )
{
    return valid( key ) ? items[key].getvalue() : T();
}

template < class T , template< class > class C >
bool ItemManager<T,C>::valid( std::string key )
{
    if( items.find(key) == items.end() ) {
        std::cout << " -- IM Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        return false;
    }
    return true;
}

#endif // KUCMSItemManagerHeader_hh


/**
 * @file VPCDisjunction.cpp
 * @author A. M. Kazachkov
 * @date 2019-12-26
 */
#include "VPCDisjunction.hpp"
#include <limits>
#include <cmath> // abs

/****************** PUBLIC  **********************/
/** Handle parameters */
void VPCDisjunction::setParams(const VPCParameters& param) {
  this->params = param;
} /* setParams */

/** Param constructor */
VPCDisjunction::VPCDisjunction(const VPCParameters& param) {
  initialize(NULL, &param);
} /* param constructor */

/** Copy and param constructor */
VPCDisjunction::VPCDisjunction(const VPCDisjunction& source, const VPCParameters& param) {
  initialize(&source, &param);
} /* copy & param constructor */

/** Default constructor */
VPCDisjunction::VPCDisjunction() {
  initialize(NULL, NULL);
} /* default constructor */

/** Copy constructor */
VPCDisjunction::VPCDisjunction(const VPCDisjunction& source) {
  initialize(&source, NULL);
} /* copy constructor */

/** Destructor */
VPCDisjunction::~VPCDisjunction() {
} /* destructor */

/** Assignment operator */
VPCDisjunction& VPCDisjunction::operator=(const VPCDisjunction& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Set up the disjunction class as new (except the timer pointer, and do not reset params) */
void VPCDisjunction::setupAsNew() {
  Disjunction::setupAsNew();
} /* setupAsNew */

/****************** PROTECTED **********************/
void VPCDisjunction::initialize(const VPCDisjunction* const source,
    const VPCParameters* const params) {
  Disjunction::initialize(source);
  if (params != NULL) {
    setParams(*params);
  }
  if (source != NULL) {
    if (params == NULL) {
      setParams(source->params);
    }
    this->timer = source->timer;
  } else {
    this->timer = NULL; // do *not* add this to setupAsNew
  }
} /* initialize */

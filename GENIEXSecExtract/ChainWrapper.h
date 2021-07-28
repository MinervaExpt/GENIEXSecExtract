#ifndef CHAINWRAPPER_H
#define CHAINWRAPPER_H

#include "TreeWrapper.h"

class ChainWrapper : public TreeWrapper
{
public:

  // Default ROOT constructor. Don't use
  ChainWrapper() {};

  explicit ChainWrapper(const char* name);

  virtual ~ChainWrapper() { delete tree; }

  int Add(const char* glob);
};
#endif

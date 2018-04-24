// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_DPG_RADIATIVE_TRANSFER_PASSKEY_HH
#define DUNE_DPG_RADIATIVE_TRANSFER_PASSKEY_HH

namespace Dune {
  template <typename T>
  class PassKey {
    friend T;

    PassKey() {}
    PassKey(const PassKey&) = default;
    PassKey& operator=(const PassKey&) = delete;
  };
}

#endif

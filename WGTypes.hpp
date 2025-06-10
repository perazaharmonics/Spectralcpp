 /* * Description:
 * * This file contains the ccommon waveguide objects used in the wg 
 * * pgsynth project.
 * *
 * * Author:
 * * JEP J. Enrique Peraza, P&G Labs, LLC
 * *
 */
#pragma once
#include <cstddef>
#include <array>
#include <cstdint>

namespace dsp::wg
{

using Sample=float;

// Forward references for visitor pattern
struct DelayBranch;
struct ScatteringJunction;

// A generic processing node
struct Node 
{
    virtual void Propagate(size_t n) noexcept = 0;
    virtual ~Node(void)=default;
};

/* Edge-table entry so InstrumentGraph knows connections */
struct Conn {Node* src; Node* dst;};

} // namespace dsp::wg
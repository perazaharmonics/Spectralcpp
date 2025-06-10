 /* * Description:
 * * This file contains the ccommon waveguide objects used in the wg 
 * * pgsynth project.
 * *
 * * Author:
 * * JEP J. Enrique Peraza
 * *
 */
#pragma once
#include <cstddef>
#include <array>
#include <cstdint>

namespace dsp::wg
{
template<typename T>
using Sample=T;

// Forward references for visitor pattern
template<size_t MaxLen>
struct DelayBranch;
template<size_t N>
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

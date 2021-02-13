#include "fragile.h"
#include <assert.h>


fragile::H::HH* fragile::H::HH::copy() {
    return (HH*)0;
}

fragile::I fragile::gI;

void fragile::fglobal(int, double, char) {
    /* empty; only used for doc-string testing */
}

namespace fragile {

    class Kderived : public K {
    public:
        virtual ~Kderived();
    };

} // namespace fragile

fragile::Kderived::~Kderived() {}

fragile::K::~K() {}

fragile::K* fragile::K::GimeK(bool derived) {
    if (!derived) return this;
    else {
        static Kderived kd;
        return &kd;
    }
};

fragile::K* fragile::K::GimeL() {
    static L l;
    return &l;
}

fragile::L::~L() {}


int fragile::create_handle(OpaqueHandle_t* handle) {
    *handle = (OpaqueHandle_t)0x01;
    return 0x01;
}

int fragile::destroy_handle(OpaqueHandle_t handle, intptr_t addr) {
    if ((intptr_t)handle == addr)
        return 1;
    return 0;
}


// for signal -> exception testing
void fragile::segfault() {
    int* i = 0; *i = 42;
}

void fragile::sigabort() {
    assert(0);
}


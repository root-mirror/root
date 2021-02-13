#include "doc_helper.h"

void DocHelper::throw_an_error(int i) {
    if (i) throw SomeError{"this is an error"};
    throw SomeOtherError{"this is another error"};
}

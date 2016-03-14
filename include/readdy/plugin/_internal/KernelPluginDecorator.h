//
// Created by clonker on 14.03.16.
//

#ifndef READDY2_MAIN_KERNELPLUGINDECORATOR_H
#define READDY2_MAIN_KERNELPLUGINDECORATOR_H


#include <readdy/plugin/Kernel.h>
#include <boost/dll/shared_library.hpp>

namespace readdy {
    namespace plugin {
        namespace _internal {
            class KernelPluginDecorator : public readdy::plugin::Kernel {
            protected:
                readdy::plugin::Kernel reference;
                boost::dll::shared_library lib;

            public:
                KernelPluginDecorator(readdy::plugin::Kernel &&reference, boost::dll::shared_library &&lib);
                ~KernelPluginDecorator() {
                    BOOST_LOG_TRIVIAL(debug) << "destroying decorator of "<< getName();
                    // TODO: weird.
                    reference.~Kernel();
                    lib.unload();
                }

                virtual const std::string getName() const override;
            };
        }
    }
}


#endif //READDY2_MAIN_KERNELPLUGINDECORATOR_H
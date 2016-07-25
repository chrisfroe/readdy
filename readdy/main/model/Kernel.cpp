/**
 * << detailed description >>
 *
 * @file Kernel.cpp.cpp
 * @brief << brief description >>
 * @author clonker
 * @date 02.05.16
 */

#include <readdy/common/make_unique.h>
#include <readdy/model/Kernel.h>
#include <atomic>

namespace readdy {
    namespace model {

        struct Kernel::Impl {
            /**
             * The name of the kernel.
             */
            std::string name;
            /**
             * todo
             */
            std::unique_ptr<signal_t> signal = std::make_unique<signal_t>();
            /**
             * todo
             */
            std::unique_ptr<_internal::ObservableFactory> observableFactory;
            /**
             * todo
             */
            std::unordered_map<ObservableBase *, boost::signals2::shared_connection_block> observableBlocks{};
        };

        const std::string &Kernel::getName() const {
            return pimpl->name;
        }

        Kernel::Kernel(const std::string &name) : pimpl(std::make_unique<Kernel::Impl>()) {
            pimpl->name = name;
            pimpl->observableFactory = std::make_unique<_internal::ObservableFactory>(this);
            std::srand((unsigned int) std::time(0));
        }

        Kernel::~Kernel() {
        }

        boost::signals2::scoped_connection Kernel::connectObservable(ObservableBase *const observable) {
            boost::signals2::scoped_connection connection (pimpl->signal.get()->connect(std::bind(&ObservableBase::callback, observable, std::placeholders::_1)));
            boost::signals2::shared_connection_block block{connection, false};
            pimpl->observableBlocks[observable] = block;
            return connection;
        }

        std::tuple<std::unique_ptr<readdy::model::ObservableWrapper>, boost::signals2::scoped_connection> Kernel::registerObservable(const ObservableType &observable, unsigned int stride) {
            auto&& wrap = std::make_unique<ObservableWrapper>(this, observable, stride);
            auto&& connection = connectObservable(wrap.get());
            return std::make_tuple(std::move(wrap), std::move(connection));
        }

        readdy::model::_internal::ObservableFactory &Kernel::getObservableFactory() const {
            return *pimpl->observableFactory;
        }

        void Kernel::evaluateObservables(readdy::model::time_step_type t) {
            for (auto &&e : pimpl->observableBlocks) {
                if (e.first->getStride() > 0 && t % e.first->getStride() != 0) {
                    e.second.block();
                } else {
                    e.second.unblock();
                }
            }
            (*pimpl->signal)(t);
        }

        void Kernel::evaluateAllObservables(readdy::model::time_step_type t) {
            for (auto &&e : pimpl->observableBlocks) {
                e.second.unblock();
            }
            (*pimpl->signal)(t);
        }

        void Kernel::deconnectObservable(ObservableBase *const observable) {
            pimpl->observableBlocks.erase(observable);
        }

        std::vector<std::string> Kernel::getAvailablePotentials() const {
            return std::vector<std::string>();
        }

        void Kernel::addParticle(const std::string &type, const Vec3 &pos) {
            getKernelStateModel().addParticle({pos[0], pos[1], pos[2], getKernelContext().getParticleTypeID(type)});
        }

        std::unique_ptr<readdy::model::potentials::Potential> Kernel::createPotential(std::string &name) const {
            return getPotentialFactory().createPotential(name);
        }


        Kernel &Kernel::operator=(Kernel &&rhs) = default;

        Kernel::Kernel(Kernel &&rhs) = default;
    }
}



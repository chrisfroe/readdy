/**
 * << detailed description >>
 *
 * @file Utils.h
 * @brief << brief description >>
 * @author clonker
 * @date 13.07.16
 */

#ifndef READDY_TESTING_UTILS_H
#define READDY_TESTING_UTILS_H

#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/log/trivial.hpp>
#include <readdy/common/Utils.h>

namespace readdy {
    namespace testing {

        inline std::vector<std::string> getKernelsToTest() {
            #ifdef READDY_KERNELS_TO_TEST
            std::vector<std::string> kernels = readdy::utils::split(std::string(READDY_KERNELS_TO_TEST), ',');
            #else
            std::vector<std::string> kernels {"SingleCPU"};
            #endif
            return kernels;
        }

        inline std::string getPluginsDirectory() {
            // test for several environment variables
            const std::string envs[] {"CONDA_ENV_PATH", "CONDA_PREFIX", "PREFIX"};
            const char *env = nullptr;
            for(auto&& key : envs) {
                env = std::getenv(key.c_str());
                if(env) {
                    BOOST_LOG_TRIVIAL(trace) << "Using env-variable for plugin dir prefix " << key << "=" << env;
                    break;
                }
            }
            std::string pluginDir = "lib/readdy_plugins";
            if (env) {
                auto _env = std::string(env);
                if (!boost::algorithm::ends_with(env, "/")) {
                    _env = _env.append("/");
                }
                pluginDir = _env.append(pluginDir);
            } else {
                BOOST_LOG_TRIVIAL(trace) << "no environment variables found that indicate plugins dir.";
            }
            return pluginDir;
        }
    }
}
#endif //READDY_TESTING_UTILS_H
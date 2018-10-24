pipeline {
//    parameters {string(defaultValue: '', description: 'api-token', name: 'API_TOKEN')
//                string(defaultValue: '', description: 'Token for readthedocs', name: 'RTD_TOKEN')
//	            string(defaultValue: '', description: 'buildable phid', name: 'TARGET_PHID')
//	            string(defaultValue: 'docker_debian_testing', description: 'docker file to use', name: 'DOCKERFILE')
//                string(defaultValue: '', description: 'Commit id', name: 'COMMIT_ID')
//    } 

    agent any

//    environment {
//        OMPI_MCA_plm = 'isolated'
//        OMPI_MCA_btl = 'tcp,self'
//    }  
    options {
        disableConcurrentBuilds()
    }
    stages {
        stage ('wipe build') {
//            when {
//                anyOf{
//                    changeset glob: "**/*.cmake"
//                    changeset glob: "**/CMakeLists.txt"
//                }
//            } 
            steps {
                sh ' rm -rf build_*'
            }
        }
        stage ('configure') {
            parallel {
                stage ('testing') {
                    agent any
		    steps {
		        sh "apt-get -qq update && apt-get -y -qq install libeigen3-dev libboost-test-dev python3-dev python3-numpy python3-scipy python3-breathe python3-sphinx python3-pytest git doxygen g++ clang cmake"
                        configure('testing')
                    }
                }
            }
        }
        stage ('build') {
            parallel {
                stage ('testing') {
                    agent any
		    steps {
		        sh "apt-get -qq update && apt-get -y -qq install libeigen3-dev libboost-test-dev python3-dev python3-numpy python3-scipy python3-breathe python3-sphinx python3-pytest git doxygen g++ clang cmake"
                        build('testing')
                    }
                }
            }
        }

        stage ('Warnings gcc') {
            steps {
                warnings(consoleParsers: [[parserName: 'GNU Make + GNU C Compiler (gcc)']])
            }
        }

        stage ('test') {
            parallel {
                ////////////////////////////////////////////////////
                stage ('testing') {
                    agent any
		    steps {
		        sh "apt-get -qq update && apt-get -y -qq install libeigen3-dev libboost-test-dev python3-dev python3-numpy python3-scipy python3-breathe python3-sphinx python3-pytest git doxygen g++ clang cmake"
                        run_test('testing', 'g++')
                    }
                    post {
                        always {
                            collect_test_results('testing', 'g++')
                        }
                    }
                }
            }
        }
    }


//    post {
//        always {
//            createartifact()
//        }
//	    success {
//            //send_fail_pass('pass')
//            //trigger_readthedocs()
//        }
//
//	    failure {
//            //send_fail_pass('fail')
//        }
//    }
}


def configure(container_name) {
    def BUILD_DIR = "build_${container_name}"
    for (CXX_COMPILER in ["g++", "clang++"]) {
        sh """
mkdir -p ${BUILD_DIR}_${CXX_COMPILER}
cd ${BUILD_DIR}_${CXX_COMPILER}
CXX=${CXX_COMPILER} cmake -DCMAKE_BUILD_TYPE:STRING=Release -DRUNNING_IN_CI=ON ..
"""
    }
}

def build(container_name) {
    def BUILD_DIR = "build_${container_name}"
    for (CXX_COMPILER in ["g++", "clang++"]) {
        sh "make -C ${BUILD_DIR}_${CXX_COMPILER}"
    }

}

def run_test(container_name, cxx_compiler) {
    def BUILD_DIR = "build_${container_name}"
    sh "cd ${BUILD_DIR}_${cxx_compiler} && ctest || true"
}

def send_fail_pass(state) {
    sh """
set +x
curl https://c4science.ch/api/harbormaster.sendmessage \
-d api.token=${API_TOKEN} \
-d buildTargetPHID=${TARGET_PHID} \
-d type=${state}
"""
}

def trigger_readthedocs() {
    sh """
set +x
curl -X POST \
-d "token=${RTD_TOKEN}" \
-d "branches=master" \
https://readthedocs.org/api/v2/webhook/muspectre/26537/

curl -X POST \
-d "token=${RTD_TOKEN}" \
https://readthedocs.org/api/v2/webhook/muspectre/26537/
"""
}


def collect_test_results(container_name, cxx_compiler) {
    def BUILD_DIR = "build_${container_name}"
    junit "${BUILD_DIR}_${cxx_compiler}/test_results*.xml"
}

def createartifact() {
    sh """ set +x
curl https://c4science.ch/api/harbormaster.createartifact \
-d api.token=${API_TOKEN} \
-d buildTargetPHID=${TARGET_PHID} \
-d artifactKey="Jenkins URI" \
-d artifactType=uri \
-d artifactData[uri]=${BUILD_URL} \
-d artifactData[name]="View Jenkins result" \
-d artifactData[ui.external]=1
"""

}


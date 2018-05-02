#ifndef DWC_H
#define DWC_H

class DWC {
   private:
       static constexpr double m_C[36] = {
               0.9545870548474235, 0.0, 0.0, 0.2979321310596127, 0.0, 0.0, 
               0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 
               0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 
               0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 
               0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
               0.2979321310596127, 0.0, 0.0, -0.9545870548474235, 0.0, 0.0
             };

   public:
           DWC() {
               C = Eigen::Map<const Eigen::MatrixXd>(m_C,6,6);
           }
           virtual ~DWC() {}
           Eigen::MatrixXd C;
};

#endif
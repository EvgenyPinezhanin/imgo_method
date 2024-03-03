#ifndef _I_CHARACTERISTIC_OPT_METHOD_H_
#define _I_CHARACTERISTIC_OPT_METHOD_H_

namespace opt {
    template<typename TrialType>
    class ICharacteristicOptMethod {
    protected:
        size_t t;

        virtual void insertInSorted(const TrialType &trial) = 0;
        virtual void calcCharacteristic() = 0;

    public:
        ICharacteristicOptMethod() : t(0) {};
    };
}

#endif // _I_CHARACTERISTIC_OPT_METHOD_H_

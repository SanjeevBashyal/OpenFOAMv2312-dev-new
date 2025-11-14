/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    contactModel

Description
    methods to solve particle-particle contact

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/

#ifndef subVolume_H
#define subVolume_H

#include "treeBoundBox.H"
#include "volumeType.H"

#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct ibSubVolumeInfo
{
    ibSubVolumeInfo():
        volumeType_(volumeType::UNKNOWN)
    {}

    ibSubVolumeInfo(volumeType vType)
    {
        switch(vType)
        {
            case volumeType::UNKNOWN:
            case volumeType::MIXED:
                volumeType_ = volumeType::UNKNOWN;
                break;
            case volumeType::INSIDE:
                volumeType_ = volumeType::INSIDE;
                break;
            case volumeType::OUTSIDE:
                volumeType_ = volumeType::OUTSIDE;
                break;
        }
    }

    // Copy constructor
    ibSubVolumeInfo(const ibSubVolumeInfo& other):
        volumeType_(other.volumeType_),
        shapesIn_(other.shapesIn_.valid() ? new labelList(*other.shapesIn_) : nullptr)
    {}

    // Assignment operator
    ibSubVolumeInfo& operator=(const ibSubVolumeInfo& other)
    {
        if (this != &other)
        {
            volumeType_ = other.volumeType_;
            shapesIn_.reset(other.shapesIn_.valid() ? new labelList(*other.shapesIn_) : nullptr);
        }
        return *this;
    }

    volumeType volumeType_;
    autoPtr<labelList> shapesIn_;
};

/*---------------------------------------------------------------------------*\
                        Class subVolume Declaration
\*---------------------------------------------------------------------------*/
class subVolume : public treeBoundBox
{
    private:

    std::shared_ptr<subVolume> parentSV_;
    List<subVolume> childSubVolumes_;

    ibSubVolumeInfo cVolumeInfo_;
    ibSubVolumeInfo tVolumeInfo_;

    bool isEdge_;

    public:

    // Constructors ------------------------------------------------
    subVolume();

    subVolume
    (
        const boundBox bb
    );

    subVolume
    (
        const boundBox bb,
        const std::shared_ptr<subVolume> parentSV,
        const volumeType cVolumeType,
        const volumeType tVolumeType
    );

    // Destructors -------------------------------------------------
    ~subVolume();

    // Copy constructor
    subVolume(const subVolume& other);

    // Assignment operator
    subVolume& operator=(const subVolume& other);

    // Member Functions --------------------------------------------

    ibSubVolumeInfo& cVolumeInfo();
    ibSubVolumeInfo& tVolumeInfo();

    List<subVolume>& childSubVolumes();
    bool hasChildSubVolumes() const;

    std::shared_ptr<subVolume>& parentSV();

    ibSubVolumeInfo& getVolumeInfo(bool cIb);
    const ibSubVolumeInfo& getVolumeInfo(bool cIb) const;

    void setAsEdge();
    bool isEdge() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

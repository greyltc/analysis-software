from datetime import datetime
from pydantic.dataclasses import dataclass
from dataclasses import field
from typing import List, Tuple, Optional, Dict
from pydantic import StrictStr


class MyConfig:
    validate_assignment = True  # runtime data type checking
    arbitrary_types_allowed = True


@dataclass(config=MyConfig)  # type: ignore
class Intensity:
    value: float = field(
        metadata={"info": "Intensity value for this cell",
                  "units": "suns"})
    raw_diodes: Optional[Tuple[int, int]] = field(
        default=None,
        metadata={"info": "Raw intensity diode ADC counts",
                  "units": "ADC Counts"})
    values: Optional[Tuple[float, float]] = field(
        default=None,
        metadata={"info": "Number of suns from each diode",
                  "units": "(suns, suns)"})
    diode_cal: Optional[Tuple[int, int]] = field(
        default=None,
        metadata={"info": "Raw intensity diode ADC counts for one sun",
                  "units": "ADC Counts per sun"})


@dataclass(config=MyConfig)  # type: ignore
class MeasurementData:
    measurements: List[Tuple[float, float, float, int]] = field(
        metadata={"info": "List of measurement data points",
                  "units": "[(volts, amps, seconds, status bytes)]"})

    @property
    def v(self):
        """list of voltage values"""
        col = 0
        return [x[col] for x in self.measurements]

    @property
    def i(self):
        """list of current values"""
        col = 1
        return [x[col] for x in self.measurements]

    @property
    def t(self):
        """list of time values"""
        col = 2
        return [x[col] for x in self.measurements]

    @property
    def s(self):
        """status list"""
        col = 3
        return [x[col] for x in self.measurements]

    @property
    def comp(self):
        """list of True/False if the sourcemeter was in compliance"""
        comp_bitmask = (1 << 14)
        stats = self.s
        return [bool(x & comp_bitmask) for x in stats]

    @property
    def ovp(self):
        """list of True/False if the sourcemeter was in overvoltage"""
        comp_bitmask = (1 << 13)
        stats = self.s
        return [bool(x & comp_bitmask) for x in stats]

    @property
    def v_valid(self):
        """list of voltage values from valid measurements"""
        v = self.v
        a = self.ovp
        b = self.comp
        return [x for x, y in zip(v, a and b) if y]

    @property
    def i_valid(self):
        """list of current values from valid measurements"""
        i = self.i
        a = self.ovp
        b = self.comp
        return [x for x, y in zip(i, a and b) if y]

    @property
    def t_valid(self):
        """list of time values from valid measurements"""
        t = self.t
        a = self.ovp
        b = self.comp
        return [x for x, y in zip(t, a and b) if y]


@dataclass(config=MyConfig)  # type: ignore
class Curve:
    roi_name: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Dictionary key for the corresponding ROI"})
    forward: Optional[bool] = field(
        default=None,
        metadata={"info": "True if the sweep direction was I_sc --> V_oc"})
    v_oc: Optional[float] = field(
        default=None,
        metadata={"info": "open circuit voltage",
                  "units": "volt"})
    i_sc: Optional[float] = field(
        default=None,
        metadata={"info": "short circuit current",
                  "units": "amp"})
    p_max: Optional[float] = field(
        default=None,
        metadata={"info": "maximum power",
                  "units": "watt"})
    v_mpp: Optional[float] = field(
        default=None,
        metadata={"info": "max power point voltage",
                  "units": "volt"})
    ff: Optional[float] = field(
        default=None,
        metadata={"info": "fill factor",
                  "units": "percent"})
    fit_function: Optional[]
    # fit: Optional[fit] = field(
    #    default=None,
    #    metadata={"info": "Fit class",)
    # TODO: Figure out the best way to represent an arbitrary fit function
    # and all its fit parameters and their confidence intervals


@dataclass(config=MyConfig)  # type: ignore
class Cell:
    operator: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Person who ran the experiment"})
    time_stamp: Optional[datetime] = field(
        default=None,
        metadata={"info": "datetime object for measurement start time"})
    pcb_firmware_revision: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Firmware version in measurement PCB"})
    control_software_revision: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Version of data collection software"})
    file_format_revision: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Version of the data storage format"})
    run_description: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Measurement description"})
    analog_io: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Analog measurement source"})
    intensity: Optional[Intensity] = field(
        default=None,
        metadata={"info": "Intesity class"})
    substrate_UID: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "unique identifier for substrate"})
    adapter_board_r: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "adapter baord resistor value",
                  "units": "Ohm"})
    sample_layout: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "sample layout type"})
    user_vars: Optional[List[Dict[StrictStr, str]]] = field(
        default=None,
        metadata={"info": "{experimental variable: value} dict list",
                  "units": "suns"})
    data: Optional[MeasurementData] = field(
        default=None,
        metadata={"info": "Measurement data class"})
    v_oc: Optional[float] = field(
        default=None,
        metadata={"info": "steady state open circuit voltage",
                  "units": "volt"})
    i_sc: Optional[float] = field(
        default=None,
        metadata={"info": "steady state short circuit current",
                  "units": "amp"})
    p_max: Optional[float] = field(
        default=None,
        metadata={"info": "steady state maximum power",
                  "units": "watt"})
    v_mpp: Optional[float] = field(
        default=None,
        metadata={"info": "max power point voltage",
                  "units": "volt"})
    ff: Optional[float] = field(
        default=None,
        metadata={"info": "fill factor",
                  "units": "percent"})
    area: Optional[float] = field(
        default=None,
        metadata={"info": "illuminated device area",
                  "units": "m^2"})
    dark_area: Optional[float] = field(
        default=None,
        metadata={"info": "electrical device area",
                  "units": "m^2"})
    rois: Optional[List[Dict[StrictStr, Tuple[int, int]]]] = field(
        default=None,
        metadata={"info": "{roi name: (start index, length)} dict list"})
    substrate_position: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Measurement setup substrate slot"})
    pixel_position: Optional[StrictStr] = field(
        default=None,
        metadata={"info": "Pixel position on substrate"})
    ivs: Optional[List[Curve]] = field(
        default=None,
        metadata={"info": "List of I-V curve classes"})


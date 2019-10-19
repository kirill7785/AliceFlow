object FormResidualSpallart_Allmares: TFormResidualSpallart_Allmares
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'residual'
  ClientHeight = 385
  ClientWidth = 705
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object Chart1: TChart
    Left = 0
    Top = 0
    Width = 705
    Height = 385
    Title.Text.Strings = (
      'Residual')
    LeftAxis.Automatic = False
    LeftAxis.AutomaticMaximum = False
    LeftAxis.AutomaticMinimum = False
    LeftAxis.Logarithmic = True
    LeftAxis.Maximum = 1.000000000000000000
    LeftAxis.Minimum = 0.000000000000010000
    View3D = False
    TabOrder = 0
    DefaultCanvas = 'TGDIPlusCanvas'
    ColorPaletteIndex = 13
    object Series1: TFastLineSeries
      Legend.Text = 'X-Vel'
      LegendTitle = 'X-Vel'
      LinePen.Color = 10708548
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {00010000000E71AC8B11BA7940}
    end
    object Series2: TFastLineSeries
      Legend.Text = 'Y-Vel'
      LegendTitle = 'Y-Vel'
      LinePen.Color = 3513587
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {00010000002B7E8CB9B9F06D40}
    end
    object Series3: TFastLineSeries
      Legend.Text = 'Z-Vel'
      LegendTitle = 'Z-Vel'
      LinePen.Color = 1330417
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {0001000000DC5ED298F6127540}
    end
    object Series4: TFastLineSeries
      Legend.Text = 'continity'
      LegendTitle = 'continity'
      LinePen.Color = 11048782
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {0001000000D061BE7C64B48040}
    end
    object Series5: TFastLineSeries
      Legend.Text = 'nut'
      LegendTitle = 'nut'
      LinePen.Color = 7028779
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {0001000000984682AFB9117640}
    end
  end
  object Timer1: TTimer
    OnTimer = Timer1Timer
    Left = 632
    Top = 320
  end
end

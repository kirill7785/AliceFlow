object FormResidualSATemp: TFormResidualSATemp
  Left = 0
  Top = 0
  Caption = 'Residual'
  ClientHeight = 435
  ClientWidth = 699
  Color = clMoneyGreen
  CustomTitleBar.CaptionAlignment = taCenter
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnResize = FormResize
  PixelsPerInch = 96
  TextHeight = 13
  object Chart1: TChart
    Left = 0
    Top = 0
    Width = 697
    Height = 433
    Title.Text.Strings = (
      'TChart')
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
      Title = 'X-Vel'
      LinePen.Color = 10708548
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {0001000000827E3AF93AD88E40}
    end
    object Series2: TFastLineSeries
      Title = 'Y-Vel'
      LinePen.Color = 3513587
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {00010000001013B94A6F1E8E40}
    end
    object Series3: TFastLineSeries
      Title = 'Z-Vel'
      LinePen.Color = 1330417
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {00010000003D56F64B886D8F40}
    end
    object Series4: TFastLineSeries
      Title = 'continity'
      LinePen.Color = 11048782
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {0001000000062AC27922588840}
    end
    object Series5: TFastLineSeries
      Title = 'Temperature'
      LinePen.Color = 7028779
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {0001000000801CC7DCD62F8D40}
    end
    object Series6: TFastLineSeries
      Title = 'nut'
      LinePen.Color = 6519581
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
      Data = {00010000003D58C3B28D6A8A40}
    end
  end
  object TimerSATemp: TTimer
    OnTimer = TimerSATempTimer
    Left = 648
    Top = 360
  end
end

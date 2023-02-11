object FormResidual_Lagtry_Menter_Temp: TFormResidual_Lagtry_Menter_Temp
  Left = 0
  Top = 0
  Caption = 'Residual'
  ClientHeight = 464
  ClientWidth = 723
  Color = clBtnFace
  CustomTitleBar.CaptionAlignment = taCenter
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  OnResize = FormResize
  PixelsPerInch = 96
  TextHeight = 13
  object Chart1: TChart
    Left = 0
    Top = 0
    Width = 721
    Height = 465
    Title.Text.Strings = (
      'Residual')
    LeftAxis.Automatic = False
    LeftAxis.AutomaticMaximum = False
    LeftAxis.AutomaticMinimum = False
    LeftAxis.Logarithmic = True
    LeftAxis.Maximum = 1438.985600000001000000
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
    end
    object Series2: TFastLineSeries
      Legend.Text = 'Y-Vel'
      LegendTitle = 'Y-Vel'
      LinePen.Color = 3513587
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object Series3: TFastLineSeries
      Legend.Text = 'Z-Vel'
      LegendTitle = 'Z-Vel'
      LinePen.Color = 1330417
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object Series4: TFastLineSeries
      Legend.Text = 'continity'
      LegendTitle = 'continity'
      LinePen.Color = 11048782
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object Series5: TFastLineSeries
      Legend.Text = 'Temperature'
      LegendTitle = 'Temperature'
      LinePen.Color = 7028779
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object Series6: TFastLineSeries
      Legend.Text = 'k'
      LegendTitle = 'k'
      LinePen.Color = 6519581
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object Series7: TFastLineSeries
      Legend.Text = 'omega'
      LegendTitle = 'omega'
      LinePen.Color = 919731
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object Series8: TFastLineSeries
      Legend.Text = 'gamma'
      LegendTitle = 'gamma'
      LinePen.Color = 6144242
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object Series9: TFastLineSeries
      Legend.Text = 'ReTheta'
      LegendTitle = 'ReTheta'
      LinePen.Color = 10401629
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
  end
  object Timer1: TTimer
    OnTimer = Timer1Timer
    Left = 624
    Top = 352
  end
end

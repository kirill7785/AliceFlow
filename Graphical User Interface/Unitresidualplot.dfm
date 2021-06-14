object Formresidual: TFormresidual
  Left = 192
  Top = 124
  Caption = 'residual'
  ClientHeight = 561
  ClientWidth = 769
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
  object cht1: TChart
    Left = 0
    Top = 0
    Width = 769
    Height = 561
    BackWall.Brush.Style = bsClear
    Title.Text.Strings = (
      'Residual')
    LeftAxis.Automatic = False
    LeftAxis.AutomaticMaximum = False
    LeftAxis.AutomaticMinimum = False
    LeftAxis.AxisValuesFormat = '#0.###E-0'
    LeftAxis.ExactDateTime = False
    LeftAxis.Logarithmic = True
    LeftAxis.Maximum = 10.000000000000000000
    LeftAxis.Minimum = 0.000000000000000001
    View3D = False
    TabOrder = 0
    DefaultCanvas = 'TGDIPlusCanvas'
    ColorPaletteIndex = 13
    object lnsrsSeries1: TLineSeries
      SeriesColor = clRed
      Title = 'X-Vel'
      Brush.BackColor = clDefault
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object lnsrsSeries2: TLineSeries
      SeriesColor = clGreen
      Title = 'Y-Vel'
      Brush.BackColor = clDefault
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object lnsrsSeries3: TLineSeries
      SeriesColor = clBlue
      Title = 'Z-Vel'
      Brush.BackColor = clDefault
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
    object lnsrsSeries4: TLineSeries
      SeriesColor = clOlive
      Title = 'Continity'
      Brush.BackColor = clDefault
      Pointer.InflateMargins = True
      Pointer.Style = psRectangle
      XValues.Name = 'X'
      XValues.Order = loAscending
      YValues.Name = 'Y'
      YValues.Order = loNone
    end
  end
  object Timer1: TTimer
    OnTimer = Timer1Timer
    Left = 704
    Top = 384
  end
  object ApplicationEvents1: TApplicationEvents
    OnMessage = ApplicationEvents1Message
    Left = 32
    Top = 16
  end
end

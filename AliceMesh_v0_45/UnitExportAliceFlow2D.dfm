object FormExportAliceFlow2D: TFormExportAliceFlow2D
  Left = 0
  Top = 0
  Caption = 'Export Alice Flow 2D'
  ClientHeight = 246
  ClientWidth = 418
  Color = clMoneyGreen
  CustomTitleBar.CaptionAlignment = taCenter
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object GroupBox1: TGroupBox
    Left = 8
    Top = 8
    Width = 281
    Height = 153
    Caption = 'Vertical vibrations'
    TabOrder = 0
    object Panel1: TPanel
      Left = 16
      Top = 48
      Width = 257
      Height = 89
      TabOrder = 0
      object Label1: TLabel
        Left = 16
        Top = 16
        Width = 46
        Height = 13
        Caption = 'amplitude'
      end
      object Label2: TLabel
        Left = 208
        Top = 16
        Width = 31
        Height = 13
        Caption = 'micron'
      end
      object Label3: TLabel
        Left = 16
        Top = 56
        Width = 49
        Height = 13
        Caption = 'frequency'
      end
      object Label4: TLabel
        Left = 216
        Top = 48
        Width = 12
        Height = 13
        Caption = 'Hz'
      end
      object Editamplitude: TEdit
        Left = 68
        Top = 13
        Width = 121
        Height = 21
        TabOrder = 0
        Text = '100'
      end
      object Editfrequency: TEdit
        Left = 71
        Top = 48
        Width = 121
        Height = 21
        TabOrder = 1
        Text = '20'
      end
    end
    object CheckBox1: TCheckBox
      Left = 16
      Top = 24
      Width = 97
      Height = 17
      Caption = 'active'
      Checked = True
      State = cbChecked
      TabOrder = 1
      OnClick = CheckBox1Click
    end
  end
  object ButtonRun: TButton
    Left = 320
    Top = 200
    Width = 75
    Height = 25
    Caption = 'Run'
    TabOrder = 1
    OnClick = ButtonRunClick
  end
end

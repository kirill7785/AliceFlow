object FormSpeedInitialization: TFormSpeedInitialization
  Left = 0
  Top = 0
  AutoSize = True
  Caption = 'Speed Initialization'
  ClientHeight = 281
  ClientWidth = 417
  Color = clMoneyGreen
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object PanelInitializationMain: TPanel
    Left = 0
    Top = 0
    Width = 417
    Height = 281
    Color = clMoneyGreen
    ParentBackground = False
    TabOrder = 0
    object Label1: TLabel
      Left = 16
      Top = 19
      Width = 109
      Height = 13
      Caption = 'Velocity X - component'
    end
    object Label2: TLabel
      Left = 16
      Top = 59
      Width = 109
      Height = 13
      Caption = 'Velocity Y - component'
    end
    object Label3: TLabel
      Left = 16
      Top = 104
      Width = 109
      Height = 13
      Caption = 'Velocity Z - component'
    end
    object Label4: TLabel
      Left = 320
      Top = 24
      Width = 17
      Height = 13
      Caption = 'm/s'
    end
    object Label5: TLabel
      Left = 320
      Top = 64
      Width = 17
      Height = 13
      Caption = 'm/s'
    end
    object Label6: TLabel
      Left = 320
      Top = 104
      Width = 17
      Height = 13
      Caption = 'm/s'
    end
    object EditVx: TEdit
      Left = 176
      Top = 16
      Width = 121
      Height = 21
      TabOrder = 0
      Text = '0.0'
    end
    object EditVy: TEdit
      Left = 176
      Top = 56
      Width = 121
      Height = 21
      TabOrder = 1
      Text = '0.0'
    end
    object EditVz: TEdit
      Left = 176
      Top = 101
      Width = 121
      Height = 21
      TabOrder = 2
      Text = '0.0'
    end
    object Button1: TButton
      Left = 264
      Top = 152
      Width = 75
      Height = 25
      Caption = 'Apply'
      TabOrder = 3
      OnClick = Button1Click
    end
  end
end
